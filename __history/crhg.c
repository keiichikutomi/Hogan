/*ARCLMCR FOR WIN32 SINCE 2023.09.11.KEIICHIKUTOMI.*/
/*LAST CHANGE:2023.09.17.*/

int arclmCR(struct arclmframe *af,int idinput);
double *rotationvct(double **rmtx);
double **rotationmtx(double *rvct);
double **spinmtx(double *rvct);
double **spinfittermtx(double *eform);
double **projectionmtx(double *eform,double *edisp,double **G);
double **jacobimtx(double *edisp,double *estress,double **M);
double **assemshelltmtx(double **estiff,double *eform,double *edisp,double *estress,double **T,FILE *fmtx);
void symmetricmtx(double **estiff,int msize);
void updaterotation(double *ddisp,double *gvct,int nnode);
double *extractshelldisplacement(struct oshell shell,double *ddisp);
double *extractshelllocalcoord(struct oshell shell,double *gdisp);
void extractdeformation(double *eform,double *edisp,int nnod);

int arclmCR(struct arclmframe *af,int idinput)
{
  DWORDLONG memory0,memory1,memory2;

  FILE *fin,*fout,*fonl,*ffig,*ffig2,*fbcl,*fmtx;               /*FILE 8 BYTES*/
  double *iform,*ddisp,*ddisp2,*dreact;
  double *lastddisp;                 /*FOR DICIDING THE DIRECTION OF PREDICTOR*/

  struct memoryelem *melem;
  struct memoryshell *mshell;

  char dir[]=DIRECTORY;                                       /*DATA DIRECTORY*/
  char s[80],string[400];

  int i,ii,jj;
  int nnode,nelem,nshell,nsect,nreact;
  int nlap,laps;                                                   /*LAP COUNT*/
  long int loffset,msize,fnode,nline;
  long int time;

  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*g,*p;/*GLOBAL MATRIX*/
  double gg,kk;
  double *weight;                                /*ARC-LENGTH WEIGHT PARAMETER*/
  double *gvct,*gvct2,*evct,*pvct,*dpvct,*tarf,*truef,*ue,*up,*re,*rp;


  /*pvct:EXTERNAL FORCE VECTOR*/
  /*dpvct:INCREMENTAL EXTERNAL FORCE VECTOR*/
  /*tarf:TARGET INTERNAL FORCE VECTOR*/
  /*truef:INTERNAL FORCE VECTOR*/
  /*evct:UNBALANCED INTERNAL FORCE VECTOR*/

  double **drccos,**Ke,**Kt,**DBe;

  double determinant,sign,lastsign,safety,dsafety;
  double arcsum,arclength,arcsign;
  double *lastpivot;
  double **mode;/*normalized LDLmode*/
  double *norm,*dm;/*norm & pivot of LDLmode*/
  double lambda;
  double upue,upup;

  double *gdisp,*edisp;                         /*DISPLACEMENT OF EACH ELEMENT*/
  double *gform,*eform;                      /*INITIAL COORDINATION OF ELEMENT*/
  double *gstress,*estress;                                /*STRESS OF ELEMENT*/
  double *shellstress;                      /*σx,σy,τxy,Mx,My,Mxy OF EACH NODE*/

  double **T,**Tt;
  double **G,**H,**P;

  long int *loffset2,*m;

  /***FOR ARC-LENGTH INCREMENTAL***/
  int iteration;
  double residual;
  int bucklingflag=0;
  int maxiteration=20;
  double tolerance=0.00001;
  double eps;

  clock_t t0,t1,t2;

  struct osect *sects;
  struct onode *nodes,*ninit;
  struct owire elem;
  struct owire *elems;
  struct oshell shell;
  struct oshell *shells;
  struct oconf *confs;

  memory0=availablephysicalmemoryEx("INITIAL:");   /*MEMORY AVAILABLE*/

  fin=fgetstofopen(dir,"r",idinput);              /*OPEN INPUT FILE*/
  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/

  fonl=fopen("hognon.onl","w"); 		   /*ITELATION OUTPUT FILE*/
  ffig=fopen("hognon.fig","w");                /*FIGURE OUTPUT FILE*/
  ffig2=fopen("hognon.fig2","w");
  fbcl=fopen("hognon.bcl","w");
  fmtx=fopen("hognon.mtx","w");

  t0=clock();                                        /*CLOCK BEGIN.*/

  getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);   /*GET INCREMENT PARAM*/

  inputinitII(fin,&nnode,&nelem,&nshell,&nsect);              /*INPUT INITIAL.*/

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

  gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
		malloc(msize*sizeof(struct gcomponent));
  gvct=(double *)malloc(msize*sizeof(double));           /*GLOBAL VECTOR.*/
  pvct=(double *)malloc(msize*sizeof(double));           /*BASE FORCE VECTOR.*/
  evct=(double *)malloc(msize*sizeof(double));           /*UNBALANCED FORCE VECTOR.*/
  dpvct=(double *)malloc(msize*sizeof(double));          /*INCREMENT FORCE VECTOR.*/
  tarf=(double *)malloc(msize*sizeof(double));           /*TARGET FORCE VECTOR.*/
  truef=(double *)malloc(msize*sizeof(double));          /*TRUE FORCE VECTOR.*/
  ue=(double *)malloc(msize*sizeof(double));
  up=(double *)malloc(msize*sizeof(double));
  weight=(double *)malloc((msize+1)*sizeof(double));          /*ARC-LENGTH WEIGHT.*/
  lastpivot=(double *)malloc(msize*sizeof(double));       /*PIVOT SIGN OF TANGENTIAL STIFFNESS.*/


  if(gmtx==NULL || gvct==NULL || pvct==NULL || evct==NULL || dpvct==NULL || tarf==NULL || truef==NULL || up==NULL || ue==NULL || weight==NULL) return 0;
  for(i=0;i<msize;i++)
  {
	(gmtx+i)->down=NULL;
	*(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
	*(evct+i)=0.0;
	*(pvct+i)=0.0;
	*(dpvct+i)=0.0;
	*(tarf+i)=0.0;
	*(truef+i)=0.0;
	*(ue+i)=0.0;
	*(up+i)=0.0;
	*(weight+i)=0.0;
  }


  free(af->sects);
  free(af->nodes);
  free(af->ninit);
  free(af->elems);
  free(af->shells);
  free(af->confs);
  free(af->ddisp);
  free(af->melem);
  free(af->mshell);

  sects=(struct osect *)malloc(nsect*sizeof(struct osect));
  if(sects==NULL) return 0;
  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(nodes==NULL) return 0;
  ninit=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(ninit==NULL) return 0;
  elems=(struct owire *)malloc(nelem*sizeof(struct owire));
  if(elems==NULL) return 0;
  shells=(struct oshell *)malloc(nshell*sizeof(struct oshell));
  if(shells==NULL) return 0;
  confs=(struct oconf *)malloc(msize*sizeof(struct oconf));
  if(confs==NULL) return 0;
  ddisp=(double *)malloc(6*nnode*sizeof(double));
  if(ddisp==NULL) return 0;
  lastddisp=(double *)malloc(6*nnode*sizeof(double));
  if(lastddisp==NULL) return 0;
  iform=(double *)malloc(6*nnode*sizeof(double));
  if(iform==NULL) return 0;
  
  melem=(struct memoryelem *)
		malloc(nelem*sizeof(struct memoryelem));
  if(melem==NULL) return 0;
  mshell=(struct memoryshell *)
		malloc(nshell*sizeof(struct memoryshell));
  if(mshell==NULL) return 0;


  af->sects=sects;
  af->nodes=nodes;
  af->ninit=ninit;
  af->elems=elems;
  af->shells=shells;
  af->confs=confs;
  af->ddisp=ddisp;                     /*DISPLACEMENT:6 DIRECTIONS.*/
  af->melem=melem;
  af->mshell=mshell;

  inputtexttomemory(fin,af);        /*READY TO READ LONG REACTIONS.*/
  nnode=af->nnode;
  nelem=af->nelem;
  nshell=af->nshell;
  nsect=af->nsect;
  nreact=af->nreact;

  initialform(nodes,ddisp,nnode);           /*ASSEMBLAGE FORMATION.*/
  initialform(nodes,iform,nnode);              /*INITIAL FORMATION.*/

  initialelem(elems,melem,nelem);            /*ASSEMBLAGE ELEMENTS.*/
  initialshell(shells,mshell,nshell);         /*ASSEMBLAGE ELEMENTS.*/

  dreact=(double *)malloc(nreact*sizeof(double));       /*REACTION.*/
  af->dreact=dreact;
  initialreact(fin,dreact,nreact);     /*ASSEMBLAGE LONG REACTIONS.*/

  assemconfII(confs,pvct,nnode);

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/
  /*drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,0,0,255);*/                      /*DRAW GLOBAL AXIS.*/

  nlap=1;
  iteration=1;
  residual=0.0;

  /*ARC-LENGTH METHOD SETTING.*/
  fnode=321;
  for(i=0;i<nnode;i++)
  {
	if((nodes+i)->code==fnode)
	{
	  *(weight+6*i+2)=1.0;
	}
  }
  *(weight+msize)=1.0;

  while(nlap<=laps)
  {
	fprintf(fonl,"LAP:%5ld/%5ld ITERATION:%5ld\n",nlap,laps,iteration);
	af->nlaps=nlap;
	arclength=dsafety;

	for(i=1;i<=msize;i++)/*MATRIX INITIALIZATION.*/
	{
	  g=(gmtx+(i-1))->down;   /*NEXT OF DIAGONAL.*/
	  while(g!=NULL) 	      /*CLEAR ROW.*/
	  {
		p=g;
		g=g->down;
		free(p);
	  }

	  ginit.m=(unsigned short int)i;
	  *(gmtx+(i-1))=ginit;
	  *(truef+(i-1))=0.0;			 /*GLOBAL VECTOR INITIALIZATION.*/
	}
	comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

	if(iteration==1)
	{
	  fprintf(fout,"\"DISPLACEMENT\"\n");
	  outputdisp(ddisp,fout,nnode,nodes);                    /*FORMATION OUTPUT.*/
	  fprintf(fout,"\"STRESS\"\n");
	}

	for(i=1;i<=nshell;i++)
	{
	  inputshell(shells,mshell,i-1,&shell);
	  shell.sect=(shells+i-1)->sect;                      /*READ SECTION DATA.*/

	  for(ii=0;ii<shell.nnod;ii++)
	  {
		inputnode(iform,shell.node[ii]);
	  }
	  gform=extractshelldisplacement(shell,iform);                      /*{Xg}*/
	  eform=extractshelllocalcoord(shell,gform);                        /*{Xe}*/

	  for(ii=0;ii<shell.nnod;ii++)
	  {
		inputnode(ddisp,shell.node[ii]);
	  }
	  gdisp=extractshelldisplacement(shell,ddisp);                   /*{Xg+Ug}*/
	  edisp=extractshelllocalcoord(shell,gdisp);                     /*{Xe+Ue}*/


	  extractdeformation(eform,edisp,shell.nnod);                       /*{Ue}*/

	  drccos=shelldrccos(shell);                                      /*DRCCOS*/
	  T=transmatrixIII(drccos,shell.nnod);         /*TRANSFORMATION MATRIX[T].*/

	  if(iteration==1)
	  {
		DBe=(double **)malloc(18*sizeof(double *));
		for(ii=0;ii<18;ii++)
		{
		  *(DBe+ii)=(double *)malloc(18*sizeof(double));
		  for(jj=0;jj<18;jj++)
		  {
			*(*(DBe+ii)+jj)=0.0;                                          /*INITIAL.*/
		  }
		}
	  }
	  else
	  {
		DBe=NULL;
	  }

	  Ke=assemshellemtx(shell,drccos,DBe);       /*ELASTIC MATRIX OF SHELL[K].*/

	  estress=matrixvector(Ke,edisp,6*shell.nnod);              /*{Fe}=[K]{Ue}*/
										 /*estress:ELEMENT INTERNAL FORCE{Fe}.*/

	  if(iteration==1)
	  {
		shellstress=matrixvector(DBe,edisp,6*shell.nnod);
		for(ii=0;ii<shell.nnod;ii++)                            /*UPDATE STRESS.*/
		{
		  for(jj=0;jj<6;jj++)
		  {
			shell.stress[ii][jj]=*(shellstress+6*ii+jj);
		  }
		}
		outputshellstress(shell,shellstress,fout);
		free(shellstress);
		freematrix(DBe,18);
	  }

	  Kt=assemshelltmtx(Ke,eform,edisp,estress,T,fmtx);
													  /*TANGENTIAL MATRIX[Kt].*/
										  /*PROJECTION of estress[Pt][Ht]{Fe}.*/

	  Kt=transformationIII(Kt,T,6*shell.nnod);

	  symmetricmtx(Kt,6*shell.nnod);
										  /*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
	  assemgstiffnessII(gmtx,Kt,&shell);
									 /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/
	  Tt=matrixtranspose(T,6*shell.nnod);                              /*[Tt].*/

	  gstress=matrixvector(Tt,estress,6*shell.nnod);
										  /*gstress:GLOBAL INTERNAL FORCE{Fg}.*/

	  for(ii=0;ii<shell.nnod;ii++)
	  {
		for(jj=0;jj<6;jj++)
		{
		  *(truef+6*(shell.node[ii]->loff)+jj)+=*(gstress+6*ii+jj);
		}
	  }

	  freematrix(drccos,3);
	  freematrix(T,18);
	  freematrix(Tt,18);
	  freematrix(Ke,18);
	  freematrix(Kt,18);

	  free(estress);
	  free(gstress);
	  free(eform);
	  free(gform);
	  free(edisp);
	  free(gdisp);

	  if(iteration==1 && (wdraw.childs+1)->hdcC!=NULL)   /*DRAW DEFORMED ELEMENT.*/
	  {
		drawglobalshell((wdraw.childs+1)->hdcC,
						(wdraw.childs+1)->vparam,
						*af,shell,255,255,255,
								  255,255,255,0,ONSCREEN/*,i*/);
	  }

	}

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
	  }
	  /*fprintf(fonl,"\"REACTION\"\n");*/
	  /*outputreaction(gmtx,gvct,nodes,confs,dreact,fonl,nnode);*/
	}


	overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/

	/*PREDICTOR*/
	if(iteration==1)
	{
	  bucklingflag=0;
	  residual=0;
	  for(i=0;i<msize;i++)
	  {
		*(up+i)=*(pvct+i);
		if((confs+i)->iconf==1) *(truef+ii)=0.0;
		*(evct+i)=*(tarf+i)-*(truef+i);/*evct:UNBALANCED FORCE -{E}.*/
		residual+=*(evct+i)**(evct+i);
		*(ue+i)=*(evct+i);
	  }
	  fprintf(fonl,"RESIDUAL=%12.9f\n",residual);

	  nline=croutludecomposition_arclength(gmtx,
								 up,ue,confs,
								 6*nnode,
								 &determinant,&sign,iteration);  /*[K]{U_p}={P}*/
	  sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",determinant,sign,comps);
	  fprintf(fonl,"%s\n",string);

	  if(sign<0.0)
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


		fclose(fin);
		fclose(fonl);
		fclose(fout);
		fclose(ffig);
		fclose(ffig2);

		gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/
		free(gvct);
		free(evct);
		free(pvct);
		free(dpvct);
		free(truef);
		free(tarf);
		free(ue);
		free(up);
		free(weight);

		return 1;
	  }

	  if(sign>0 && nlap!=1)/*最初の修正子計算時に前のlapからのpivotの符号変化が検出された場合座屈点精算(PIN-POINTING)開始*/
	  {
		bucklingflag=sign;
		fprintf(fbcl,"\nLAP:%d/%d ITERATION:%d\n",nlap,laps,iteration);
		m=(long int *)malloc(bucklingflag*sizeof(long int));
		jj=0;
		for(ii=0;ii<msize;ii++)
		{
			if((confs+ii)->iconf==0 && (gmtx+ii)->value<0)
			{
				*(m+jj)=ii;
				jj++;
				fprintf(fbcl,"m=%5ld\n",ii);
			}
		}
	  }

	  arcsum=*(weight+msize);
	  for(ii=0;ii<msize;ii++)
	  {
		if((confs+ii)->iconf!=1)
		{
		  if(nlap==1)
		  {
			*(weight+ii)*=1.0/(arclength*arclength**(up+ii)**(up+ii));
		  }
		  arcsum+=*(weight+ii)**(up+ii)**(up+ii);/*変位スケーリングU_p*G*U_p*/
		}
	  }
	  lambda=arclength*arclength/sqrt(arcsum);/*荷重倍率dΛ=Δ/√(U_p*weight*U_p+γ)*/
	  fprintf(fonl,"ARCSUM=%12.9f ",arcsum);
	  fprintf(fonl,"LAMBDA=%12.9f ",lambda);

	  /*予測子符号決定*/
	  if(nlap==1)
	  {
		arcsign=1.0;
		for(ii=0;ii<msize;ii++)
		{
		  *(lastddisp+ii)=*(ddisp+ii);
		}
	  }
	  else
	  {
		arcsign=0.0;
		for(ii=0;ii<msize;ii++)
		{
		  arcsign+=(*(ddisp+ii)-*(lastddisp+ii))**(up+ii);
		  *(lastddisp+ii)=*(ddisp+ii);
		}
		fprintf(fonl,"ARCSIGN=%12.5f\n",arcsign);
		arcsign/=abs(arcsign);
	  }


	  for(ii=0;ii<msize;ii++)
	  {
		if((confs+ii)->iconf!=1)
		{
		  *(gvct+ii)*=*(up+ii)*lambda*arcsign;            /*gvct:増分変位dΛU_p*/
		  *(dpvct+ii)=*(pvct+ii)*lambda*arcsign;           /*dpvct:増分荷重dΛP*/
		  *(tarf+ii)+=*(dpvct+ii);                 /*tarf:増分後の目標荷重P+dP*/
		}
	  }

	  if(fonl!=NULL)
	  {
		fprintf(fonl,"\"NEW TARGET FORCE {F}\"\n");
		for(i=0;i<nnode;i++)
		{
		  fprintf(fonl,"NODE:%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",(nodes+i)->code,
				  *(tarf+(6*i+0)),*(tarf+(6*i+1)),*(tarf+(6*i+2)),
				  *(tarf+(6*i+3)),*(tarf+(6*i+4)),*(tarf+(6*i+5)));
		}
	  }
	}
	/*CORRECTOR*/
	if(iteration!=1)
	{
	  residual=0;
	  eps=1e-6;
	  for(i=0;i<msize;i++)
	  {
		*(up+i)=*(pvct+i);
		if((confs+i)->iconf==1) *(truef+i)=0.0;
		*(evct+i)=*(tarf+i)-*(truef+i);   /*unblf:UNBALANCED FORCE -{E}.*/
		residual+=*(evct+i)**(evct+i);    /*residual*/
		*(ue+i)=*(evct+i);
	  }
	  fprintf(fonl,"RESIDUAL=%12.9f\n",residual);


	  if(fonl!=NULL)
	  {
		fprintf(fonl,"\"TRUE FORCE {F-dF}:AFTER MODIFY CONF\"\n");
		for(i=0;i<nnode;i++)
		{
		  fprintf(fonl,"NODE:%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",(nodes+i)->code,*(truef+(6*i+0)),*(truef+(6*i+1)),*(truef+(6*i+2)),*(truef+(6*i+3)),*(truef+(6*i+4)),*(truef+(6*i+5)));
		}
		fprintf(fonl,"\"UNBALANCED FORCE {dF}\"\n");
		for(i=0;i<nnode;i++)
		{
		  fprintf(fonl,"NODE:%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",(nodes+i)->code,*(evct+(6*i+0)),*(evct+(6*i+1)),*(evct+(6*i+2)),*(evct+(6*i+3)),*(evct+(6*i+4)),*(evct+(6*i+5)));
		}
	  }

	  nline=croutludecomposition_arclength(gmtx,
								 up,ue,confs,
								 6*nnode,
								 &determinant,&sign,iteration);        /*[K]{U_p}={P}*/
	  sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",determinant,sign,comps);
	  fprintf(fonl,"%s\n",string);

	  if(sign<0.0)
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

		/*laptime("\0",t0);*/

		fclose(fin);
		fclose(fonl);
		fclose(fout);
		fclose(ffig);
		fclose(ffig2);

		gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/
		free(gvct);
		free(evct);
		free(pvct);
		free(dpvct);
		free(truef);
		free(tarf);
		free(ue);
		free(up);

		return 1;
	  }

	  if(sign>0 && bucklingflag==0)/*修正子計算時に前のlapからのpivotの符号変化が検出された場合座屈点精算(PIN-POINTING)開始*/
	  {
		bucklingflag=sign;
		fprintf(fbcl,"\nLAP:%d/%d ITERATION:%d\n",nlap,laps,iteration);
		m=(long int *)malloc(bucklingflag*sizeof(long int));
		jj=0;
		for(ii=0;ii<msize;ii++)
		{
			if((confs+ii)->iconf==0 && (gmtx+ii)->value<0)
			{
				*(m+jj)=ii;
				jj++;
				fprintf(fbcl,"m=%5ld\n",ii);
			}
		}
	  }


	  if(bucklingflag>0)  /*PIN-POINTING BIGIN*/
	  {
		/*FOR DIRECTIONAL DERIVATIVE OF TANGENTIAL STIFFNESS MATRIX*/
		fprintf(fbcl,"ITERATION:%d\n",iteration);
		mode=(double **)malloc(bucklingflag*sizeof(double *));
		for(ii=0;ii<bucklingflag;ii++)
		{
		  *(mode+ii)=(double *)malloc(msize*sizeof(double));
		}
		norm=(double *)malloc(bucklingflag*sizeof(double));
		dm=(double *)malloc(bucklingflag*sizeof(double));


		/*for(ii=0;ii<msize;ii++)
		{
		  fprintf(ffig2,"\nlastpivot=%12.9f",*(lastpivot+ii));
		  fprintf(ffig2,"\npivot=%12.9f",(gmtx+ii)->value);
		}*/

		LDLmode(gmtx,confs,m,bucklingflag,mode,norm,dm,msize);
		for(ii=0;ii<bucklingflag;ii++)
		{
		  fprintf(fbcl,"mode=%5ld norm=%12.9f dm=%12.9f\n",*(m+ii),*(norm+ii),*(dm+ii));
		}
		fprintf(fbcl,"\"LDL EIGEN VECTOR:\"\n");
		for(i=0;i<nnode;i++)
		{
		  fprintf(fbcl,"NODE:%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",(nodes+i)->code,
		  *(*(mode+0)+(6*i+0)),*(*(mode+0)+(6*i+1)),*(*(mode+0)+(6*i+2)),
		  *(*(mode+0)+(6*i+3)),*(*(mode+0)+(6*i+4)),*(*(mode+0)+(6*i+5)));
		}

		ddisp2=(double *)malloc(6*nnode*sizeof(double));
		gvct2=(double *)malloc(6*nnode*sizeof(double));

		if(ddisp2==NULL || gvct2==NULL) return 0;
		for(i=0;i<msize;i++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
		{
          if((confs+ii)->iconf!=1)
		  {
			*(gvct2+i)=eps**(*(mode+0)+i);
		  }
		  *(ddisp2+i)=*(ddisp+i);
		}
		updaterotation(ddisp2,gvct2,nnode);

		re=(double *)malloc(msize*sizeof(double));
		rp=(double *)malloc(msize*sizeof(double));

		for(i=1;i<=nshell;i++)
		{
		  inputshell(shells,mshell,i-1,&shell);
		  shell.sect=(shells+i-1)->sect;                      /*READ SECTION DATA.*/

		  for(ii=0;ii<shell.nnod;ii++)
		  {
			inputnode(iform,shell.node[ii]);
		  }
		  gform=extractshelldisplacement(shell,iform);                  /*{Xg}*/
		  eform=extractshelllocalcoord(shell,gform);                    /*{Xe}*/

		  for(ii=0;ii<shell.nnod;ii++)
		  {
			inputnode(ddisp2,shell.node[ii]);
		  }
		  gdisp=extractshelldisplacement(shell,ddisp2);              /*{Xg+Ug}*/
		  edisp=extractshelllocalcoord(shell,gdisp);                 /*{Xe+Ue}*/

		  extractdeformation(eform,edisp,shell.nnod);                   /*{Ue}*/

		  drccos=shelldrccos(shell);                                  /*DRCCOS*/
		  T=transmatrixIII(drccos,shell.nnod);     /*TRANSFORMATION MATRIX[T].*/

		  Ke=assemshellemtx(shell,drccos,NULL);  /*ELASTIC MATRIX OF SHELL[K].*/

		  estress=matrixvector(Ke,edisp,6*shell.nnod);          /*{Fe}=[K]{Ue}*/
										 /*estress:ELEMENT INTERNAL FORCE{Fe}.*/
		  Kt=assemshelltmtx(Ke,eform,edisp,estress,T,fmtx);
													  /*TANGENTIAL MATRIX[Kt].*/
		  Kt=transformationIII(Kt,T,6*shell.nnod);

		  symmetricmtx(Kt,6*shell.nnod);
										  /*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
		  loffset2=(long int *)malloc(18*sizeof(long int));
		  for(ii=0;ii<3;ii++)
		  {
			for(jj=0;jj<6;jj++)
			{
			  *(loffset2+(6*ii+jj))=6*(shell.node[ii]->loff)+jj;
			}
		  }
		  for(ii=0;ii<18;ii++)
		  {
			for(jj=0;jj<18;jj++)
			{
			   *(re+*(loffset2+ii))+=*(*(Kt+ii)+jj)**(ue+*(loffset2+jj));
			   *(rp+*(loffset2+ii))+=*(*(Kt+ii)+jj)**(up+*(loffset2+jj));
			}
		  }
		  freematrix(drccos,3);
		  freematrix(T,18);
		  freematrix(Ke,18);
		  freematrix(Kt,18);
		  free(loffset2);
		  free(estress);
		  free(eform);
		  free(gform);
		  free(edisp);
		  free(gdisp);
		}
		upue=*(dm+0)/(*(norm+0)**(norm+0));

		fprintf(fbcl,"PIN-POINTED EIGENVALUE=%12.9f\n",upue);

		upue*=eps;
		upup=0;
		for(i=0;i<msize;i++)
		{
		  if((confs+i)->iconf!=1)
		  {
			*(rp+i)-=*(pvct+i);
			*(re+i)-=*(evct+i);
			upup+=*(rp+i)**(*(mode+0)+i);
			upue+=*(re+i)**(*(mode+0)+i);
		  }
		}
		fprintf(fonl,"upue=%12.9f upup=%12.9f\n",upue,upup);
		lambda=-upue/upup;


		freematrix(mode,bucklingflag);
		free(norm);
		free(dm);
		free(ddisp2);
		free(gvct2);
		free(re);
		free(rp);
	  }
	  else/*修正子計算時に前のlapからのpivotの符号変化が検出されない場合、次のiterationへ行き計算続行*/
	  {

		/*不平衡残差最小法*/
		upue=0;
		upup=*(weight+msize);
		for(ii=0;ii<=msize;ii++)
		{
		  if((confs+ii)->iconf!=1)
		  {
			upue+=*(weight+ii)*(*(up+ii))*(*(ue+ii));/*U_p*G*δU_e*/
			upup+=*(weight+ii)*(*(up+ii))*(*(up+ii));/*U_p*G*U_p*/
			/*fprintf(fonl,"gvct=%12.9f ubf=%12.9f\n",*(gvct+ii),*(ubf+ii));*/
		  }
		}
		fprintf(fonl,"upue=%12.9f upup=%12.9f\n",upue,upup);
		lambda=-upue/upup;
	  }


	  fprintf(fonl,"LAMBDA=%12.9f\n",lambda);
	  for(ii=0;ii<=msize;ii++)
	  {
		if((confs+ii)->iconf!=1)
		{
		  *(gvct+ii)=*(up+ii)*lambda+*(ue+ii);/*gvct:増分変位δU_e+δΛU_p*/
		  *(dpvct+ii)=*(pvct+ii)*lambda;/*dpvct:増分荷重δΛP*/
		  *(tarf+ii)+=*(dpvct+ii);/*tarf:増分後の荷重P+dP*/
		}
	  }
	  if(fonl!=NULL)
	  {
		fprintf(fonl,"\"NEW TARGET FORCE {F}:\"\n");
		for(i=0;i<nnode;i++)
		{
		  fprintf(fonl,"NODE:%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",(nodes+i)->code,
		  *(tarf+(6*i+0)),*(tarf+(6*i+1)),*(tarf+(6*i+2)),
		  *(tarf+(6*i+3)),*(tarf+(6*i+4)),*(tarf+(6*i+5)));
		}
	  }
	}

	updaterotation(ddisp,gvct,nnode);                      /*FORMATION UPDATE.*/
	updaterotation(pvct,gvct,nnode);                      /*FORMATION UPDATE.*/

	for(ii=0;ii<nnode;ii++)
	{
	  if(iteration==1)
	  {
		if((nodes+ii)->code==fnode && ffig2!=NULL)
		{
		  fprintf(ffig2,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f {RESISUAL}= %16.12f {DETERMINANT}= %16.12f {SIGN}= %5.0f {BUCKLINGFLAG}= %1d\n",
				  nlap,laps,(nodes+ii)->code,*(truef+(6*ii+2)),*(ddisp+6*ii+2)-*(iform+6*ii+2),residual,determinant,sign,bucklingflag);
		}
	  }
	  if((nodes+ii)->code==fnode && ffig!=NULL)
	  {
		fprintf(ffig,"LAP: %3d / %3d ITERATION: %2d NODE %3d {Fz}= %18.12f {U}= %16.12f {RESISUAL}= %16.12f {DETERMINANT}= %16.12f {SIGN}= %5.0f {BUCKLINGFLAG}= %1d\n",
				  nlap,laps,iteration,(nodes+ii)->code,*(truef+(6*ii+2)),*(ddisp+6*ii+2)-*(iform+6*ii+2),residual,determinant,sign,bucklingflag);
	  }
	}

	/***KUTOMI FOR ARC-LENGTH INCREMENTAL***/
	if( (residual<tolerance || iteration==maxiteration)&& iteration!=1/*&& bucklingflag==0*/)
	{
	  sprintf(string,"LAP:%5ld/%5ld FINISHED.",nlap,laps);
	  laptime(string,t0);
	  nlap++;
	  iteration=0;
	}
	iteration++;


    while(!GetAsyncKeyState(VK_LBUTTON))  /*LEFT CLICK TO CONTINUE.*/
	{
	  if(GetAsyncKeyState(VK_RBUTTON))      /*RIGHT CLICK TO ABORT.*/
	  {
		fclose(fin);
		gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/
		free(gvct);
		free(pvct);
		free(dpvct);
		free(truef);
		free(tarf);
		free(ue);
		free(up);

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
  }                                        /*REPEAT UNTIL INSTABLE.*/

  if((wdraw.childs+1)->hdcC!=NULL && ddisp!=NULL)                 /*DRAW LAST FRAME.*/
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
	for(i=1;i<=nshell;i++)
	{
	  shell=*(shells+i-1);                     /*READ ELEMENT DATA.*/

	  for(ii=0;ii<shell.nnod;ii++)
	  {
		inputnode(ddisp,shell.node[ii]);
	  }
	  if(globaldrawflag==1)
	  {
		drawglobalshell((wdraw.childs+1)->hdcC,
					   (wdraw.childs+1)->vparam,
					   *af,shell,255,255,255,
							   255,255,255,0,ONSCREEN);
	  }
	}
	overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }

  fclose(fin);
  fclose(fout);
  fclose(ffig);
  fclose(ffig2);

  gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/

  af->eigenvec=(double **)malloc(1*sizeof(double *));
  *((af->eigenvec)+0)=gvct;

  errormessage(" ");
  errormessage("COMPLETED.");

  if(fonl!=NULL) fclose(fonl);

  return 0;
}/*arclmCR*/

double *rotationvct(double **rmtx)
{
  double c,s;
  double *rvct;
  double theta;
  rvct=(double *)malloc(3*sizeof(double));
  c=*(*(rmtx+0)+0)+*(*(rmtx+1)+1)+*(*(rmtx+2)+2)-1;                         /*2cos(alpha)*/
  s=sqrt((*(*(rmtx+2)+1)-*(*(rmtx+1)+2))*(*(*(rmtx+2)+1)-*(*(rmtx+1)+2))
		+(*(*(rmtx+0)+2)-*(*(rmtx+2)+0))*(*(*(rmtx+0)+2)-*(*(rmtx+2)+0))
		+(*(*(rmtx+1)+0)-*(*(rmtx+0)+1))*(*(*(rmtx+1)+0)-*(*(rmtx+0)+1)));/*2sin(alpha)>0*/
  if(c>0)
  {
	theta=atan(s/c);
  }
  else
  {
	theta=atan(s/c)+PI;
  }
  if(s!=0.0)
  {
	*(rvct+0)=theta*((*(*(rmtx+2)+1)-*(*(rmtx+1)+2)))/s;
	*(rvct+1)=theta*((*(*(rmtx+0)+2)-*(*(rmtx+2)+0)))/s;
	*(rvct+2)=theta*((*(*(rmtx+1)+0)-*(*(rmtx+0)+1)))/s;
  }
  else
  {
	*(rvct+0)=0.0;
	*(rvct+1)=0.0;
	*(rvct+2)=0.0;
  }
  return rvct;
}

double **rotationmtx(double *rvct)
{
  int i;
  double **rmtx;
  double *n;
  double theta;

  n=(double *)malloc(3*sizeof(double));
  rmtx=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(rmtx+i)=(double *)malloc(3*sizeof(double));
  }
  theta=sqrt(*(rvct+0)**(rvct+0)+*(rvct+1)**(rvct+1)+*(rvct+2)**(rvct+2));
  if(theta!=0)
  {
	for(i=0;i<3;i++)
	{
	  *(n+i)=*(rvct+i)/theta;
	}

	*(*(rmtx+0)+0)= cos(theta)       +(1-cos(theta))**(n+0)**(n+0);
	*(*(rmtx+0)+1)=-sin(theta)**(n+2)+(1-cos(theta))**(n+0)**(n+1);
	*(*(rmtx+0)+2)= sin(theta)**(n+1)+(1-cos(theta))**(n+0)**(n+2);
	*(*(rmtx+1)+0)= sin(theta)**(n+2)+(1-cos(theta))**(n+1)**(n+0);
	*(*(rmtx+1)+1)= cos(theta)       +(1-cos(theta))**(n+1)**(n+1);
	*(*(rmtx+1)+2)=-sin(theta)**(n+0)+(1-cos(theta))**(n+1)**(n+2);
	*(*(rmtx+2)+0)=-sin(theta)**(n+1)+(1-cos(theta))**(n+2)**(n+0);
	*(*(rmtx+2)+1)= sin(theta)**(n+0)+(1-cos(theta))**(n+2)**(n+1);
	*(*(rmtx+2)+2)= cos(theta)       +(1-cos(theta))**(n+2)**(n+2);
  }
  else
  {
	*(*(rmtx+0)+0)=1.0;
	*(*(rmtx+0)+1)=0.0;
	*(*(rmtx+0)+2)=0.0;
	*(*(rmtx+1)+0)=0.0;
	*(*(rmtx+1)+1)=1.0;
	*(*(rmtx+1)+2)=0.0;
	*(*(rmtx+2)+0)=0.0;
	*(*(rmtx+2)+1)=0.0;
	*(*(rmtx+2)+2)=1.0;
  }
  free(n);
  return rmtx;
}

double **spinmtx(double *rvct)
{
  int i;
  double **smtx;

  smtx=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(smtx+i)=(double *)malloc(3*sizeof(double));
  }

  *(*(smtx+0)+0)=0.0;
  *(*(smtx+0)+1)=-*(rvct+2);
  *(*(smtx+0)+2)=*(rvct+1);
  *(*(smtx+1)+0)=*(rvct+2);
  *(*(smtx+1)+1)=0.0;
  *(*(smtx+1)+2)=-*(rvct+0);
  *(*(smtx+2)+0)=-*(rvct+1);
  *(*(smtx+2)+1)=*(rvct+0);
  *(*(smtx+2)+2)=0.0;

  return smtx;
}

double **spinfittermtx(double *eform)
{
  int i,j;
  double A,len;
  double x1,x2,y1,y2;
  double **G;

  G=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(G+i)=(double *)malloc(18*sizeof(double));
	for(j=0;j<18;j++)
	{
	  *(*(G+i)+j)=0.0;
	}
  }

  x1=*(eform+6)-*(eform+0);
  y1=*(eform+7)-*(eform+1);
  x2=*(eform+12)-*(eform+0);
  y2=*(eform+13)-*(eform+1);
  A=0.5*(x1*y2-x2*y1);
  len=sqrt(x1*x1+y1*y1);

  *(*(G+0)+2)=0.5*(*(eform+12)-*(eform+6))/A;
  *(*(G+1)+2)=0.5*(*(eform+13)-*(eform+7))/A;
  *(*(G+2)+1)=-1/len;

  *(*(G+0)+8)=0.5*(*(eform+0)-*(eform+12))/A;
  *(*(G+1)+8)=0.5*(*(eform+1)-*(eform+13))/A;
  *(*(G+2)+7)= 1/len;

  *(*(G+0)+14)=0.5*(*(eform+6)-*(eform+0))/A;
  *(*(G+1)+14)=0.5*(*(eform+7)-*(eform+1))/A;

  return G;
}

double **projectionmtx(double *eform,double *edisp,double **G)
{
  int i,j,a,b;
  double *node;
  double **P,**Gu,**S,**SGu;

  P=(double **)malloc(18*sizeof(double *));
  for(i=0;i<18;i++)
  {
	*(P+i)=(double *)malloc(18*sizeof(double));
  }

  node=(double *)malloc(3*sizeof(double));

  Gu=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Gu+i)=(double *)malloc(3*sizeof(double));
  }

  for(a=0;a<3;a++)
  {
	for(i=0;i<3;i++)
	{
	  *(node+i)=*(eform+6*a+i)+*(edisp+6*a+i);
	}
	S=spinmtx(node);
	for(b=0;b<3;b++)
	{
	  for(i=0;i<3;i++)
	  {
		for(j=0;j<3;j++)
		{
		  *(*(Gu+i)+j)=*(*(G+i)+6*b+j);
		}
	  }
	  SGu=matrixmatrix(S,Gu,3);
	  for(i=0;i<3;i++)
	  {
		for(j=0;j<3;j++)
		{
		  *(*(P+6*a+3+i)+6*b+j)=-*(*(Gu+i)+j);
		  *(*(P+6*a+i)+6*b+3+j)=0.0;
		  if(a==b && i==j)
		  {
			*(*(P+6*a+i)+6*b+j)=2.0/3.0;
			*(*(P+6*a+3+i)+6*b+3+j)=1.0;
		  }
		  else if(i==j)
		  {
			*(*(P+6*a+i)+6*b+j)=-1.0/3.0;
			*(*(P+6*a+3+i)+6*b+3+j)=0.0;
		  }
		  else
		  {
			*(*(P+6*a+i)+6*b+j)=0.0;
			*(*(P+6*a+3+i)+6*b+3+j)=0.0;
		  }
		  *(*(P+6*a+i)+6*b+j)+=*(*(SGu+i)+j);

		}
	  }
	  freematrix(SGu,3);
	}
	freematrix(S,3);
  }
  free(node);
  freematrix(Gu,3);
  return P;
}

double **jacobimtx(double *edisp,double *estress,double **M)
{
  DWORDLONG memory0,memory1,memory2;
  int n,i,j;
  double theta,eta,mu,dot;
  double *thetavct,*mvct;
  double **thetaspin,**thetaspin2,**mspin,**mtheta,**mtheta2;
  double **H,**Ha,**Ma;


  H=(double **)malloc(18*sizeof(double *));
  for(i=0;i<18;i++)
  {
	*(H+i)=(double *)malloc(18*sizeof(double));
	for(j=0;j<18;j++)
	{
	 *(*(H+i)+j)=0.0;
	}
  }

  mtheta=(double **)malloc(3*sizeof(double *));                    /*{ma}{θat}*/
  for(i=0;i<3;i++)
  {
	*(mtheta+i)=(double *)malloc(3*sizeof(double));
  }

  thetavct=(double *)malloc(3*sizeof(double));
  mvct=(double *)malloc(3*sizeof(double));

  for(n=0;n<3;n++)
  {
	theta=0.0;
	dot=0.0;
	for(i=0;i<3;i++)
	{
	  *(thetavct+i)=*(edisp+6*n+3+i);
	  *(mvct+i)=*(estress+6*n+3+i);
	  theta+=*(thetavct+i)**(thetavct+i);
	  dot+=*(thetavct+i)**(mvct+i);
	}
	theta=sqrt(theta);
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		 *(*(mtheta+i)+j)=*(mvct+i)**(thetavct+j);
	  }
	}
	
	if(theta<PI/30.0)
	{
	  eta=1.0/12.0+pow(theta,2)/720.0+pow(theta,4)/30240.0+pow(theta,6)/1209600.0;
	  mu=1.0/360.0+pow(theta,2)/7560.0+pow(theta,4)/201600.0+pow(theta,6)/5987520.0;
	}
	else
	{
	  eta=(1.0-0.5*theta*(1.0/tan(0.5*theta)))/pow(theta,2);
	  mu=(theta*theta+4.0*cos(theta)+theta*sin(theta)-4.0)
		 /(4.0*pow(theta,4)*sin(0.5*theta)*sin(0.5*theta));
	}
	thetaspin=spinmtx(thetavct);
	thetaspin2=matrixmatrix(thetaspin,thetaspin,3);
	mspin=spinmtx(mvct);
	mtheta2=matrixmatrix(thetaspin2,mtheta,3);
	Ha=(double **)malloc(3*sizeof(double *));
	for(i=0;i<3;i++)
	{
	  *(Ha+i)=(double *)malloc(3*sizeof(double));
	}
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Ha+i)+j)=-0.5**(*(thetaspin+i)+j)+eta**(*(thetaspin2+i)+j);
		if(i==j)*(*(Ha+i)+j)+=1.0;
		*(*(mtheta2+i)+j)*=mu;
		*(*(mtheta2+i)+j)+=eta*(*(*(mtheta+j)+i)-2.0**(*(mtheta+i)+j))-0.5**(*(mspin+i)+j);
		if(i==j)*(*(mtheta2+i)+j)+=eta*dot;
	  }
	}

	Ma=matrixmatrix(mtheta2,Ha,3);

	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(H+6*n+3+i)+6*n+3+j)=*(*(Ha+i)+j);
		*(*(M+6*n+3+i)+6*n+3+j)=*(*(Ma+i)+j);
		if(i==j)*(*(H+6*n+i)+6*n+j)=1.0;
	  }
	}
	freematrix(thetaspin,3);
	freematrix(thetaspin2,3);
	freematrix(mspin,3);
	freematrix(mtheta2,3);
	freematrix(Ma,3);
	freematrix(Ha,3);
  }
  free(thetavct);
  free(mvct);
  freematrix(mtheta,3);
  return H;
}

double **assemshelltmtx(double **estiff,double *eform,double *edisp,double *estress,double **T,FILE *fmtx)
{
  DWORD memory0,memory1,memory2;
  int i,j,n;
  double **G,**P,**H;
  double **HP,**PtHt;
  double *nm,**spinnm,*pstress;
  double **Fnm,**Fn,**FnG,**GtFnt;
  double **tstiff,**gstiff1,**gstiff2,**gstiff3;


  tstiff=(double **)malloc(18*sizeof(double *));                              /*[M]*/
  for(i=0;i<18;i++)
  {
	*(tstiff+i)=(double *)malloc(18*sizeof(double));
  }

  gstiff3=(double **)malloc(18*sizeof(double *));                              /*[M]*/
  for(i=0;i<18;i++)
  {
	*(gstiff3+i)=(double *)malloc(18*sizeof(double));
	for(j=0;j<18;j++)
	{
	  *(*(gstiff3+i)+j)=0.0;
	}
  }


  G=spinfittermtx(eform);                             /*SPIN-FITTER MATRIX[G].*/
  P=projectionmtx(eform,edisp,G);                            /*PROJECTION MATRIX[P].*/
  H=jacobimtx(edisp,estress,gstiff3);              /*JACOBIAN MATRIX OF ROTATION[H].*/
  HP=matrixmatrix(H,P,18);                                            /*[H][P]*/
  PtHt=matrixtranspose(HP,18);
  pstress=matrixvector(PtHt,estress,18);           /*projected estress {Fp}*/
  for (i=0;i<18;i++)
  {
	*(estress+i)=*(pstress+i);
  }
  free(pstress);

  nm=(double *)malloc(3*sizeof(double));           /*projected estress of each node {n}&{m}*/

  Fnm=(double **)malloc(18*sizeof(double *));                          /*[Fnm]*/
  for(i=0;i<18;i++)
  {
	*(Fnm+i)=(double *)malloc(3*sizeof(double));
  }
  Fn=(double **)malloc(18*sizeof(double *));                            /*[Fn]*/
  for(i=0;i<18;i++)
  {
	*(Fn+i)=(double *)malloc(3*sizeof(double));
  }

  for(n=0;n<3;n++)
  {
	for(i=0;i<3;i++)
	{
	  *(nm+i)=*(estress+6*n+i);
	}
	spinnm=spinmtx(nm);
	for(i=0;i<3;i++)                                    
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Fnm+6*n+i)+j)=*(*(spinnm+i)+j);
		*(*(Fn+6*n+i)+j)=*(*(spinnm+i)+j);
	  }
	}
	freematrix(spinnm,3);

	for(i=0;i<3;i++)
	{
	  *(nm+i)=*(estress+6*n+3+i);
	}
	spinnm=spinmtx(nm);
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Fnm+6*n+3+i)+j)=*(*(spinnm+i)+j);
		*(*(Fn+6*n+3+i)+j)=0.0;
	  }
	}
	freematrix(spinnm,3);
  }


  estiff=transformationIII(estiff,HP,18);                  /*[Pt][Ht][K][H][P]*/

  gstiff1=matrixmatrixIII(Fnm,G,18,3,18);                           /*[Fnm][G]*/

  FnG=matrixmatrixIII(Fn,G,18,3,18);                                 /*[Fn][G]*/
  GtFnt=matrixtranspose(FnG,18);                                   /*[Gt][Fnt]*/
  gstiff2=matrixmatrix(GtFnt,P,18);                             /*[Gt][Fnt][P]*/

  gstiff3=transformationIII(gstiff3,P,18);                        /*[Pt][M][P]*/

  for(i=0;i<18;i++)
  {
	for(j=0;j<18;j++) 
	{
	  *(*(tstiff+i)+j)=*(*(estiff+i)+j)-*(*(gstiff1+i)+j)-*(*(gstiff2+i)+j)+*(*(gstiff3+i)+j);
	}
  }


  free(nm);

  freematrix(G,3);
  freematrix(P,18);
  freematrix(H,18);
  freematrix(PtHt,18);
  freematrix(HP,18);
  freematrix(Fnm,18);
  freematrix(Fn,18);

  freematrix(gstiff1,18);
  freematrix(FnG,18);
  freematrix(GtFnt,18);
  freematrix(gstiff2,18);
  freematrix(gstiff3,18);
  return tstiff;
}

void symmetricmtx(double **estiff,int msize)
{
  int i,j;
  for(i=0;i<msize;i++)
  {
	for(j=0;j<i;j++) 
	{
	  *(*(estiff+i)+j)=0.5*(*(*(estiff+i)+j)+*(*(estiff+j)+i));
	  *(*(estiff+j)+i)=*(*(estiff+i)+j);
	}		
  }
  return;
}

void updaterotation(double *ddisp,double *gvct,int nnode)
/*FORMATION UPDATE IF ROTATION IS FINITE.*/
{
  int i,j;
  long int loff;
  double *rvctR,*rvctL,*rvct;
  double **rmtxR,**rmtxL,**rmtx;
  rvctR=(double *)malloc(3*nnode*sizeof(double));
  rvctL=(double *)malloc(3*nnode*sizeof(double));
  for(i=0;i<nnode;i++)
  {
	for(j=0;j<3;j++)
	{
	  loff=6*i+j;
	  *(ddisp+loff)+=*(gvct+loff);
	}
	for(j=0;j<3;j++)
	{
	  loff=6*i+3+j;
	  *(rvctR+j)=*(ddisp+loff);
	  *(rvctL+j)=*(gvct+loff);
	}
	rmtxR=rotationmtx(rvctR);
	rmtxL=rotationmtx(rvctL);
	rmtx=matrixmatrix(rmtxL,rmtxR,3);
	rvct=rotationvct(rmtx);
	for(j=0;j<3;j++)
	{
	  loff=6*i+3+j;
	  *(ddisp+loff)=*(rvct+j);
	}
    freematrix(rmtxR,3);
	freematrix(rmtxL,3);
	freematrix(rmtx,3);
	free(rvct);
  }
  free(rvctR);
  free(rvctL);
  return;
}/*updaterotation*/

double *extractshelldisplacement(struct oshell shell,double *ddisp)
/*EXTRACT ELEMENT DEFORMATION{dU} FROM GLOBAL VECTOR.*/
{
  long int i,loffset;
  int n,nnod;
  double *d;

  nnod=shell.nnod;
  d=(double *)malloc(6*nnod*sizeof(double));
  if(d==NULL)
  {
	errormessage("EXTRACTDISPLACEMENT:MEMORY INSUFFICIENT.");
	return NULL;
  }
  for(n=0;n<nnod;n++)
  {
	for(i=0;i<6;i++)
	{
	  loffset=6*(shell.node[n]->loff)+i;
	  *(d+6*n+i)=*(ddisp+loffset);
	}
  }
  return d;
}/*extractshelldisplacement*/

double *extractshelllocalcoord(struct oshell shell,double *gdisp)
/*EXTRACT LOCAL ELEMENT DEFORMATION FROM GLOBAL.*/
/*UPDATE PSUEDO-ROTATION VECTOR*/
{
  long int i,loffset;
  int n,nnod;
  double *d,*r,*c,*td,*tr,*edisp;
  double **drccos,**trmtx,**rmtx;

  nnod=shell.nnod;

  edisp=(double *)malloc(6*nnod*sizeof(double));

  c=(double *)malloc(3*sizeof(double));
  d=(double *)malloc(3*sizeof(double));
  r=(double *)malloc(3*sizeof(double));

  drccos=shelldrccos(shell);

  for(i=0;i<3;i++)
  {
	*(c+i)=0.0;
	for(n=0;n<nnod;n++)
	{
	  *(c+i)+=*(gdisp+6*n+i)/nnod;
	}
  }
  for(n=0;n<nnod;n++)
  {
	for(i=0;i<3;i++)
	{
	  *(d+i)=*(gdisp+6*n+i)-*(c+i);
	  *(r+i)=*(gdisp+6*n+3+i);
	}

	td=matrixvector(drccos,d,3);
	/*location vector of each node in Local Frame from center*/

	rmtx=rotationmtx(r);
	/*rmtx:Ra*/
	trmtx=matrixmatrix(drccos,rmtx,3);
	tr=rotationvct(trmtx);

	for(i=0;i<3;i++)
	{
	  *(edisp+6*n+i)  =*(td+i);
	  *(edisp+6*n+3+i)=*(tr+i);
	}
	free(td);
	free(tr);
	freematrix(rmtx,3);
	freematrix(trmtx,3);
  }


  free(c);
  free(d);
  free(r);
  freematrix(drccos,3);

  return edisp;
}/*extractshelllocalcoord*/

void extractdeformation(double *eform,double *edisp,int nnod)
/*EXTRACT LOCAL ELEMENT DEFORMATION FROM GLOBAL.*/
/*UPDATE PSUEDO-ROTATION VECTOR*/
{
  int n,i;
  double *r,*rinit,*rvct;
  double **rmtx,**rh,**rt,**rtt;

  r    =(double *)malloc(3*sizeof(double));
  rinit=(double *)malloc(3*sizeof(double));

  for(n=0;n<nnod;n++)
  {
	for(i=0;i<3;i++)
	{
	  *(r+i)    =*(edisp+6*n+3+i);
	  *(rinit+i)=*(eform+6*n+3+i);
	}

	rh=rotationmtx(r);
	rt=rotationmtx(rinit);
	rtt=matrixtranspose(rt,3);
	rmtx=matrixmatrix(rh,rtt,3);
	rvct=rotationvct(rmtx);

	for(i=0;i<3;i++)
	{
	  *(edisp+6*n+i)  -=*(eform+6*n+i);
	  *(edisp+6*n+3+i) =*(rvct+i);
	}

	freematrix(rh,3);
	freematrix(rt,3);
	freematrix(rtt,3);
	freematrix(rmtx,3);
	free(rvct);
  }

  free(r);
  free(rinit);
  return;
}
