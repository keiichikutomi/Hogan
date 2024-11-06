
int srcan001(char fname[])
{
  int ii,jj,ie,is,ic,ia,in,na,n1,n2,n3;
  int ielem,isect,ni,nj;    /*ELEM, SECT & NODE ID OF EACH ELEMENT.*/
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

  int cmqcode[MAXCMQ];
  double Mo[MAXCMQ],Mm;
  double rate,rate1,rate2,ratema,ratemu,rateqa,ratequ;


  long int codelist[MAXSECT];
  double marate[MAXSECT],murate[MAXSECT];
  double qarate[MAXSECT],qurate[MAXSECT];

  double wrate1,wrate2; /*RATE OF WINDOW.*/

  double *cNcr;  /*Buckling Condensation Result*/
  double icNcr;
  int bclngcondensationflag;


  FILE *fin=NULL,*flist,*fsafe,*frate;
  FILE *fz=NULL,*fx=NULL,*fy=NULL;
  char txt[400],non[80],str[256];
  char strhoko[2][10]={"X","Y"};
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
  if(!strncmp(non,"projectname",11)) strcpy(prj,"projectname");
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
  /*
  if(globalmessageflag==1 &&
	 MessageBox(NULL,"Calculate Ultimate.","SRCan",MB_YESNO)==IDYES)
  {
	calc[PULTIMATE]=1;
  }
  else
  {
    calc[PULTIMATE]=0;
  }
  */

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
  fprintf(fout0,"ífñ éZíË \"S,RC,SRC\" í∑ä˙,íZä˙,èIã«\n");
  fprintf(fout0,"égópÉtÉ@ÉCÉã\n");

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,non,80);
  fprintf(fout0,"ì¸óÕÉfÅ[É^               =%s\n",non);
  fprintf(fsafe,"%s\n",non);

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,non,80);
  fprintf(fout0,"âîíºâ◊èdéûâêÕåãâ        =%s\n",non);
  fprintf(fsafe,"%s\n",non);

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,non,80);
  fprintf(fout0,"êÖïΩâ◊èdéûâêÕåãâ  Xï˚å¸ =%s\n",non);
  fprintf(fsafe,"%s\n",non);

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,non,80);
  fprintf(fout0,"êÖïΩâ◊èdéûâêÕåãâ  Yï˚å¸ =%s\n",non);
  fprintf(fsafe,"%s\n",non);

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,non,80);
  fprintf(fout0,"âºíËífñ                =%s\n",non);
  fprintf(fsafe,"%s\n",non);

  fprintf(fout0,"íPà ån tf(kN),tfm(kNm)\n");

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


  /*INITIAL SETTING*/
  jis=1.0;                                          /*1.1=JIS STEEL.*/

  if(!strcmp(prj,"projectname"))
  {
	fgetinitial2(fin,&nnode,&nelem,&nsect);   /*INPUT ARCLM INITIAL.*/
  }
  else
  {
	fgetinitial2(fin,&nnode,&nelem,&nsect);
	//fgetinitial(fin,&nnode,&nelem,&nsect,&E,&G);    /*INPUT FRAME3 INITIAL.*/
  }

  sects=(struct section *)malloc(nsect*sizeof(struct section));
  if(sects==NULL) return 0;
  nodes=(struct structnode *)malloc(nnode*sizeof(struct structnode));
  if(nodes==NULL) return 0;
  elems=(struct element *)malloc(nelem*sizeof(struct element));
  if(elems==NULL) return 0;

  if(!strcmp(prj,"projectname"))
  {
	inputfiletomemory2(fin,&nnode,&nelem,&nsect,sects,nodes,elems);
  }
  else
  {
	inputfiletomemory2(fin,&nnode,&nelem,&nsect,sects,nodes,elems);
	//inputfiletomemory(fin,&nnode,&nelem,&nsect,sects,nodes,elems);
  }

  comments(fout0);/*COMMENTS FOR TEST FILE.*/

  for(is=0;is<nsect;is++) /*SECTIONS INTO MEMORY.*/
  {
	if(getsectionform(flist,(sects+is)->code,(sects+is))==0)
	{
	  getsectionform(flist,(sects+is)->ocode,(sects+is));
	}
  }
  nsect=getcodelist(flist,codelist);
  /*nsect CHANGED.nsect IS SECTION NUM IN LIST FILE, NOT IN INL FILE.*/

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

  cNcr=(double *)malloc(nelem*sizeof(double));   /*GLOBAL VECTOR*/
  if(bclngcondensationflag)
  {
	definencr(&arc,&*cNcr);          /*Buckling Condensation Result*/
  }

  ie=0;
  while(1) /*CALCULATION*/
  {
	if(ie>=nelem) break;

	elem=*(elems+ie);            /*ALL ELEMENTS POINTING EACH SECT.*/

    if(slabcount==0 && ie<nelem-1) pair=*(elems+ie+1);
    if(slabcount==1)               pair=*(elems+ie-1);

	/*INITIALIZE STRESS*/
	initializestress(&(elem.head.x));/*SHORT X*/
	initializestress(&(elem.tail.x));
	initializestress(&(elem.head.y));/*SHORT Y*/
	initializestress(&(elem.tail.y));
	initializestress(&(elem.head.z));/*LONG*/
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

		  if     (sect1.stype==STYPE_S)     fprintf(fout0,"=Çr");
		  else if(sect1.stype==STYPE_RC)    fprintf(fout0,"=ÇqÇb");
		  else if(sect1.stype==STYPE_SRC)   fprintf(fout0,"=ÇrÇqÇb");
		  else if(sect1.stype==STYPE_PC)    fprintf(fout0,"=ÇoÇb");
		  else if(sect1.stype==STYPE_WOOD)  fprintf(fout0,"=ñÿ");
          else if(sect1.stype==STYPE_GLASS) fprintf(fout0,"=ÉKÉâÉX");
          else if(sect1.stype==STYPE_ACRYL) fprintf(fout0,"=ÉAÉNÉäÉã");
		  else if(sect1.stype==STYPE_ALUMI) fprintf(fout0,"=ÉAÉãÉ~");
		  else                              fprintf(fout0,"=  ");

		  if     (sect1.etype==COLUMN) fprintf(fout0,"íå ");
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

		  /*MATERIAL INDEPENDENT ON PERIOD.[kgf/cm2]*/
		  gmaterial.sE=2100000.0;
		  gmaterial.rE=2100000.0;
		  //gmaterial.cE=210000.0;
		  gmaterial.gE= 725000.0;
		  gmaterial.aE=  32000.0;                     /*Acryl=15000, Metacryl=32000*/
		  gmaterial.cE=gmaterial.rE/15.0;
		  gmaterial.alE=700000.0;                     /*ALUMI*/
		  gmaterial.sF=2400.0*jis;                    /*STEEL SN400*/
		  //gmaterial.alF=1750.0;                     /*ALUMI AS175*/
		  gmaterial.alF=4400.0;                       /*ALUMI A7178T6*/

		  gmaterial.Fc=240.0;                         /*CONCRETE Fc240*/



          translatesection(sect1,gmaterial,&As,&Ac,&Ar,&Yg,SX);
		  /*
		  fprintf(fout0,"As=%7.2f[cm2] ",As);
		  fprintf(fout0,"Ar=%7.2f[cm2] ",Ar);
		  fprintf(fout0,"Ac=%7.2f[cm2]\n",Ac);
		  */

		  fprintf(fout0,"âûóÕ       :        N");
          fprintf(fout0,"                Qxi                Qxj");
          fprintf(fout0,"                Qyi                Qyj");
          fprintf(fout0,"                 Mt");
          fprintf(fout0,"                Mxi                Mxj");
          fprintf(fout0,"                Myi                Myj\n");

          fprintf(fout0,"âîíºéûZ    : %8.3f(%8.2f)",
				  elem.head.z.N   ,si*elem.head.z.N);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
				  elem.head.z.Q[0],si*elem.head.z.Q[0],
                  elem.tail.z.Q[0],si*elem.tail.z.Q[0]);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.z.Q[1],si*elem.head.z.Q[1],
                  elem.tail.z.Q[1],si*elem.tail.z.Q[1]);
		  fprintf(fout0," %8.3f(%8.2f)",
				  elem.head.z.Mt  ,si*elem.head.z.Mt);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.z.M[0],si*elem.head.z.M[0],
				  elem.tail.z.M[0],si*elem.tail.z.M[0]);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)\n",
                  elem.head.z.M[1],si*elem.head.z.M[1],
				  elem.tail.z.M[1],si*elem.tail.z.M[1]);

          fprintf(fout0,"êÖïΩéûX    : %8.3f(%8.2f)",
				  elem.head.x.N   ,si*elem.head.x.N);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
				  elem.head.x.Q[0],si*elem.head.x.Q[0],
                  elem.tail.x.Q[0],si*elem.tail.x.Q[0]);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
				  elem.head.x.Q[1],si*elem.head.x.Q[1],
                  elem.tail.x.Q[1],si*elem.tail.x.Q[1]);
		  fprintf(fout0," %8.3f(%8.2f)",
				  elem.head.x.Mt  ,si*elem.head.x.Mt);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.x.M[0],si*elem.head.x.M[0],
				  elem.tail.x.M[0],si*elem.tail.x.M[0]);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)\n",
				  elem.head.x.M[1],si*elem.head.x.M[1],
				  elem.tail.x.M[1],si*elem.tail.x.M[1]);

		  fprintf(fout0,"êÖïΩéûY    : %8.3f(%8.2f)",
				  elem.head.y.N   ,si*elem.head.y.N);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
				  elem.head.y.Q[0],si*elem.head.y.Q[0],
				  elem.tail.y.Q[0],si*elem.tail.y.Q[0]);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
				  elem.head.y.Q[1],si*elem.head.y.Q[1],
				  elem.tail.y.Q[1],si*elem.tail.y.Q[1]);
		  fprintf(fout0," %8.3f(%8.2f)",
				  elem.head.y.Mt  ,si*elem.head.y.Mt);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
				  elem.head.y.M[0],si*elem.head.y.M[0],
				  elem.tail.y.M[0],si*elem.tail.y.M[0]);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)\n",
				  elem.head.y.M[1],si*elem.head.y.M[1],
				  elem.tail.y.M[1],si*elem.tail.y.M[1]);
		  fprintf(fout0,"\n");



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

		  fprintf(fout0,"í∑ä˙       : %8.3f(%8.2f)",
				  elem.head.z.N   ,si*elem.head.z.N);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
				  elem.head.z.Q[0],si*elem.head.z.Q[0],
				  elem.tail.z.Q[0],si*elem.tail.z.Q[0]);
		  fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
				  elem.head.z.Q[1],si*elem.head.z.Q[1],
				  elem.tail.z.Q[1],si*elem.tail.z.Q[1]);
		  fprintf(fout0," %8.3f(%8.2f)",
				  elem.head.z.Mt  ,si*elem.head.z.Mt);
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


		  if((sect1.stype==STYPE_S || sect1.stype==STYPE_GLASS || sect1.stype==STYPE_ACRYL || sect1.stype==STYPE_ALUMI)
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



			strcat(strQ[2][0],str);
			if(Qcrx!=0.0) sprintf(str,"           %8.3f",pstress[0]->Q[0]/Qcrx);
			else          sprintf(str,"   10.000");

			strcat(strQ[2][0],str);
			if(Qcry!=0.0) sprintf(str,"           %8.3f",pstress[0]->Q[1]/Qcry);
			else          sprintf(str,"   10.000");

			strcat(strQ[2][1],str);
			if(Mcrx!=0.0) sprintf(str,"           %8.3f",pstress[0]->M[0]/Mcrx);
			else          sprintf(str,"   10.000");

			strcat(strM[2][0],str);
			if(Mcry!=0.0) sprintf(str,"           %8.3f",pstress[0]->M[1]/Mcry);
			else          sprintf(str,"   10.000");

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



			if(!strcmp(prj,"projctname")) rate=0.0; /*PASS LONG CALCULATION.*/

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

			strcat(strQ[2][0],str);
			if(Qcry!=0.0) sprintf(str,"           %8.3f",pstress[1]->Q[1]/Qcry);
			else          sprintf(str,"   10.000");

			strcat(strQ[2][1],str);
			if(Mcrx!=0.0) sprintf(str,"           %8.3f",pstress[1]->M[0]/Mcrx);
			else          sprintf(str,"   10.000");

			strcat(strM[2][0],str);
			if(Mcry!=0.0) sprintf(str,"           %8.3f",pstress[1]->M[1]/Mcry);
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
		  else if(sect1.stype==STYPE_WOOD && sect1.sform.type==DABO) /*WOOD DABO*/
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
			  }

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



			  sprintf(str,"           %8.3f           %8.3f",
					  rate1,rate2);
			  strcat(strM[2][ia],str);

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
		  }

		  fprintf(fout0,"\n");


          /*SHORT..................................................*/
          /*MATERIAL FOR SHORT*/
          gmaterial.sft=2400.0*jis;                   /*STEEL SN400*/
          gmaterial.sfc=-gmaterial.sft;
          gmaterial.sfb=2400.0*jis; /*FOR SRC*/         /*b:BENDING*/
          gmaterial.sfs=2400.0/sqrt(3.0)*jis;             /*s:SHEAR*/
		  gmaterial.rft=3500.0*jis;           /*REINFORCEMENT SD345*/
          gmaterial.wft=3000.0;               /*REINFORCEMENT SD295*/

		  gmaterial.rfc=-gmaterial.rft;
          gmaterial.cfc=-160.0;                          /*CONCRETE*/
          gmaterial.cfs=11.1;

          for(hoko=0;hoko<=1;hoko++) /*LOAD DIRECTION X,Y.*/
		  {
			if(sect1.stype==STYPE_PC) break;  /*PC ROUTE 3a : SKIP SHORT.*/
			if(!strcmp(prj,"projectname")) break;  /*SKIP SHORT.*/
			if(!strcmp(prj,"projectname"))
			{
			  nfact=1.0;
			  mfact=1.0;
			  qfact=1.0;
			}


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
				  rate=allowablestressofflatbaralumi(PSHORT,sect1,gmaterial.alE,
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
											   h,ia,(-sN));
				  }


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



                    if(Qa>0.0) rate=fabs(sQ[in][1-ia]/Qa);
                    else       rate=10.0;

                    if(rate>rateqa) rateqa=rate;
                    if(rate>qarate[soffset-1])
                    {
                      qarate[soffset-1]=rate;
                    }
                    sprintf(str,"           %8.3f",rate);
                    strcat(strQ[2][1-ia],str);

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
          }


		  /*ULTIMATE...............................................*/

		  if(calc[PULTIMATE]==1 || sect1.stype==STYPE_PC)
		  {
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
					  Mu=ultimatebendingofsrc(elem,ia,(-uN),
											  &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
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

					  if(Qu<0.1) Qu=0.1;

					  rate=fabs(uQ[in][1-ia]/Qu);
					  if(rate>ratequ) ratequ=rate;
					  if(rate>qurate[soffset-1])
					  {
						qurate[soffset-1]=rate;
					  }
					  sprintf(str,"           %8.3f",rate);
					  strcat(strQ[2][1-ia],str);



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

		  if(frate!=NULL)
		  {
			fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f %.5f %.5f OCODE: %5d\n",
					ielem,isect,rateqa,ratequ,ratema,ratemu,sect1.ocode); /****SRCANMAX****/
		  }

		  if(ratema==10.0) dabocount++;

		}

		if(soffset && sect1.etype==WALL)/*...........SHEAR OF WALL.*/
		{
		  Qa=0.0;
		  Qu=0.0;


		  fprintf(fout0,"-------------------------------------");
		  fprintf(fout0,"-------------------------------------");
		  fprintf(fout0,"-------------------------------------");
		  fprintf(fout0,"-------------------------------------");
		  fprintf(fout0,"-------------------------------------");
		  fprintf(fout0,"-----------------\n");
		  fprintf(fout0,"ïîçﬁ:%4d ",ielem);
		  fprintf(fout0,"éní[:%3d èIí[:%3d ",ni,nj);



		  /*for original section */
		  if(sect1.code==sect1.ocode)
		  {
			fprintf(fout0,"ífñ :Å¸%3d/Å@%3d",sect1.ocode,isect);
		  }
		  else
		  {
			 fprintf(fout0,"ífñ :Å@%3d/Å¸%3d",sect1.ocode,isect);
		  }

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



            rate=fabs((fabs(lQ)+fabs(wafact*sQ[0][0]))/Qa);
            if(rate>rateqa) rateqa=rate;
            if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;
            fprintf(fout0," Qs/Qa=%7.5f",rate);

            fprintf(fout0,"\n");
          }

          if(rateqa>1.0 || ratequ>1.0)
          {
            fprintf(fsafe,"NG ELEM:%4d SECT:%4d ",ielem,isect);
            fprintf(fsafe,"MAX:Q/Qa=%.5f Q/Qu=%.5f\n",rateqa,ratequ);
          }

          if(frate!=NULL)
		  {
			fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f 0.0 0.0 OCODE: %5d\n",
					ielem,isect,rateqa,ratequ,sect1.ocode); /****SRCANMAX****/
          }
        }

        if(soffset && sect1.etype==BRACE)            /*STEEL BRACE.*/
		{
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


          if(frate!=NULL)
          {
			fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f 0.0 0.0 OCODE: %5d\n",
					ielem,isect,rateqa,ratequ,sect1.ocode); /****SRCANMAX****/
          }
        }

        if(soffset && sect1.etype==SLAB)/*...........SHEAR OF SLAB.*/
        {
          if(slabcount==0)      slabcount=1;
          else if(slabcount==1) slabcount=0;

          Qa=0.0;
		  Qu=0.0;

          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-----------------\n");
          fprintf(fout0,"ïîçﬁ:%4d ",ielem);
          fprintf(fout0,"éní[:%3d èIí[:%3d ",ni,nj);



		  /* 150429 fukushima for original section */
		  if(sect1.code==sect1.ocode)
		  {
			fprintf(fout0,"ífñ :Å¸%3d/Å@%3d",sect1.ocode,isect);
		  }
		  else
		  {
			fprintf(fout0,"ífñ :Å@%3d/Å¸%3d",sect1.ocode,isect);
		  }

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


          if(frate!=NULL)
		  {
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

  fclose(fin);  /*.inl*/
  fclose(flist);/*.lst*/
  fclose(fx);   /*.otl*/
  fclose(fy);
  fclose(fz);
  fclose(fout0);/*.tst*/
  fclose(fsafe);/*.rlt*/
  fclose(frate);/*.rat*/

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
