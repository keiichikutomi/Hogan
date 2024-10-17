


  updatestress(melem,fout,edisp,estress,ee,elem,
			   func,NULL);                               /*{f}+{df}*/

  ee:elastic
  estress : estiff(include pmtx)*edisp






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


void updateshellstress(struct memoryshell *mshell,FILE *fout,
					   double *edisp,double *dstress,double **estiff,
					   struct oshell *shell,double func[],FILE *ftxt)
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

  /*LAST LAP*/
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

  /*
  if(wsurf.hwnd!=NULL)
  {
	fseek(arc.fsurface,0L,SEEK_END);
	for(i=0;i<2;i++)
	{
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

	  if(f[i]>=pow(RADIUS,EXPONENT))
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
  */
  return;
}/*updatestress*/
