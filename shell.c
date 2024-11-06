
/*FUNCTIONS FOR SHELL*/
double*** elasticCshell(struct oshell shell);
double** shelldrccos(struct oshell shell);
double shellarea(struct oshell shell);
double* extractshelldisplacement(struct oshell shell, double* ddisp);

/*STRAIN-DISPLACEMENT MATRIX*/
double*** assemshellshape(struct oshell shell, double** drccos);
/*ELASTIC STIFFNESS MATRIX*/
double** assemshellemtx(struct oshell shell);
/*ELASTO-PLASTIC STIFFNESS MATRIX*/
double** assemshellpmtx(struct oshell shell, double*** C, double ***B);
/*MASS MATRIX*/
double** assemshellmmtx(struct oshell shell, double **drccos);
/*ELEMENT LOAD*/
double* assemshellpvct(struct oshell shell, double **drccos);
/*VOLUME*/
double shellvolume(struct oshell shell, double** drccos);
void assemshellvolume(struct oshell* shells, int nshell, double* ddisp, double* volume);


/*INPUT & OUTPUT SHELL DATA*/
void inputshell(struct oshell *shells,struct memoryshell *mshell,int offset,struct oshell *shell);
void outputshell(struct oshell *shells,int offset,struct oshell *shell);
void outputmemoryshell(struct oshell *shells,struct memoryshell *mshell,int offset);

/*FOR PLASTIC (ILYUSHIN'S STRESS RESULTANT)*/
double* ilyushin(struct oshell* shell, int ii, double* lambda, double** W);/*YIELD FUNCTION*/
double* returnmapilyushin(struct oshell* shell, int ii, double** W);/*RETURN-MAPPING*/
double** consistentCilyushin(struct oshell* shell, int ii, double** H, double Halpha, double** g, double galpha, int nstress);/*UPDATE STIFFNESS MATRIX*/

void assemshellestrain(struct oshell* shell, double*** B, double* edisp);
void assemshellestress(struct oshell* shell, double*** C);
double* assemshelleinternal(struct oshell* shell, double*** B);

double yieldstress(struct osect* sect, double alpha, double* dy, double* ddy);

extern double* extractlocalcoord(double* gform, double** drccos, double nnod);
extern double* extractshelldisplacement(struct oshell shell, double* ddisp);
extern double* extractdeformation(double* eforminit, double* eform, int nnod);

extern double** assemtmtxCR(double** Ke, double* eform, double* edisp, double* estress, double* gstress, double** T, double** HPT, int nnod);
extern void symmetricmtx(double** estiff, int msize);

extern void dbgstr(const char* str);
extern void dbgvct(double* vct, int size, int linesize, const char* str);
extern void dbgmtx(double** mtx, int rows, int cols, const char* str);


double*** elasticCshell(struct oshell shell)
{
	int i,j,ii;
	int nstress = shell.nstress;
	double*** C;
	double E = shell.sect->E;
	double t = shell.sect->area;
	double poi = shell.sect->poi;

	C = (double***)malloc(shell.ngp * sizeof(double**));
	for (ii = 0; ii < shell.ngp; ii++)
	{
		*(C+ii) = (double**)malloc(nstress * sizeof(double*));
		for (i = 0; i < nstress; i++)
		{
			*(*(C+ii) + i) = (double*)malloc(nstress * sizeof(double));
			for (j = 0; j < nstress; j++)
			{
				*(*(*(C+ii) + i) + j) = 0.0;
			}
		}
		*(*(*(C+ii)+0)+0)=E*t/(1-poi*poi);
		*(*(*(C+ii)+0)+1)=poi**(*(*(C+ii)+0)+0);
		*(*(*(C+ii)+1)+0)=*(*(*(C+ii)+0)+1);
		*(*(*(C+ii)+1)+1)=*(*(*(C+ii)+0)+0);
		*(*(*(C+ii)+2)+2)=0.5*E*t/(1+poi);

		*(*(*(C+ii)+3)+3)=E*pow(t,3)/(12.0*(1-poi*poi));
		*(*(*(C+ii)+3)+4)=poi**(*(*(C+ii)+3)+3);
		*(*(*(C+ii)+4)+3)=*(*(*(C+ii)+3)+4);
		*(*(*(C+ii)+4)+4)=*(*(*(C+ii)+3)+3);
		*(*(*(C+ii)+5)+5)=0.5*E*pow(t,3)/(12.0*(1+poi));

		if(nstress==7)*(*(*(C+ii)+6)+6)=E*t*27.0/50.0;
	}
	return C;
}


double **shelldrccos(struct oshell shell)
/*RETURN:FILM DIRECTION COSINE.*/
{
  char str[256];
  int i;
  double dl,Xx,Yx,Zx,Xy,Yy,Zy,Xz,Yz,Zz;
  double **drccos;
  struct onode n1,n2,n3;

  n1=*(shell.node[0]);
  n2=*(shell.node[1]);
  n3=*(shell.node[2]);

  drccos=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(drccos+i)=(double *)malloc(3*sizeof(double));
  }

  Xx=n2.d[GX]-n1.d[GX]; /*AXIS X*/
  Yx=n2.d[GY]-n1.d[GY];
  Zx=n2.d[GZ]-n1.d[GZ];
  dl=sqrt(Xx*Xx+Yx*Yx+Zx*Zx);
  if(dl==0.0)
  {
	sprintf(str,"N1=%ld N2=%ld N3=%ld",n1.code,n2.code,n3.code);
	MessageBox(NULL,str,"shelldrccos",MB_OK);
  }
  *(*(drccos+GX)+0)=Xx/dl;
  *(*(drccos+GX)+1)=Yx/dl;
  *(*(drccos+GX)+2)=Zx/dl;

  Xz=Yx*(n3.d[GZ]-n1.d[GZ])-Zx*(n3.d[GY]-n1.d[GY]); /*AXIS Z*/
  Yz=Zx*(n3.d[GX]-n1.d[GX])-Xx*(n3.d[GZ]-n1.d[GZ]);
  Zz=Xx*(n3.d[GY]-n1.d[GY])-Yx*(n3.d[GX]-n1.d[GX]);
  dl=sqrt(Xz*Xz+Yz*Yz+Zz*Zz);
  if(dl==0.0)
  {
	sprintf(str,"N1=%ld N2=%ld N3=%ld",n1.code,n2.code,n3.code);
	MessageBox(NULL,str,"shelldrccos",MB_OK);
  }
  *(*(drccos+GZ)+0)=Xz/dl;
  *(*(drccos+GZ)+1)=Yz/dl;
  *(*(drccos+GZ)+2)=Zz/dl;

  Xy=Yz*Zx-Zz*Yx; /*AXIS Y*/
  Yy=Zz*Xx-Xz*Zx;
  Zy=Xz*Yx-Yz*Xx;
  dl=sqrt(Xy*Xy+Yy*Yy+Zy*Zy);
  if(dl==0.0)
  {
	sprintf(str,"N1=%ld N2=%ld N3=%ld",n1.code,n2.code,n3.code);
	MessageBox(NULL,str,"shelldrccos",MB_OK);
  }
  *(*(drccos+GY)+0)=Xy/dl;
  *(*(drccos+GY)+1)=Yy/dl;
  *(*(drccos+GY)+2)=Zy/dl;

  return drccos;
}/*shelldrccos*/

double shellarea(struct oshell shell)
{
  double area,X,Y,Z,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3;
  struct onode n1,n2,n3;

  /*shell.node[i]->d[j] == *(shell.node[i]).d[j]*/
  n1=*(shell.node[0]);
  n2=*(shell.node[1]);
  n3=*(shell.node[2]);

  X1=n2.d[GX]-n1.d[GX];
  Y1=n2.d[GY]-n1.d[GY];
  Z1=n2.d[GZ]-n1.d[GZ];

  X2=n3.d[GX]-n1.d[GX];
  Y2=n3.d[GY]-n1.d[GY];
  Z2=n3.d[GZ]-n1.d[GZ];


  X=Y1*Z2-Z1*Y2;
  Y=Z1*X2-X1*Z2;
  Z=X1*Y2-Y1*X2;

  area=0.5*sqrt(X*X+Y*Y+Z*Z);

  return area;
}

double* extractshelldisplacement(struct oshell shell, double* ddisp)
/*EXTRACT ELEMENT DEFORMATION{dU} FROM GLOBAL VECTOR.*/
{
	long int i, loffset;
	int n, nnod;
	double* d;

	nnod = shell.nnod;
	d = (double*)malloc(6 * nnod * sizeof(double));
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 6; i++)
		{
			loffset = 6 * (shell.node[n]->loff) + i;
			*(d + 6 * n + i) = *(ddisp + loffset);
		}
	}
	return d;
}/*extractshelldisplacement*/

double*** assemshellshape(struct oshell shell, double** drccos)
{
  int i,j,k,ii;
  int nstress = shell.nstress;
  double *a,*b,*c;
  double *Lx,*Ly;
  double **exy;
  double area;
  double** L;
  double *aa,*bb,*cc,*dd,*ee;
  double len;
  double** Bp,** Bb,** Bt;
  double*** B;


  a=(double *)malloc(3*sizeof(double));
  b=(double *)malloc(3*sizeof(double));
  c=(double *)malloc(3*sizeof(double));
  Lx=(double *)malloc(3*sizeof(double));
  Ly=(double *)malloc(3*sizeof(double));
  exy=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)/*NODE 0-0, 0-1, 0-2*/
  {
	*(exy+i)=(double *)malloc(2*sizeof(double));
  }
  L=(double **)malloc(shell.ngp*sizeof(double *));/*TRIANGULAR COORDINATION OF INTEGRATION POINTS*/
  for(i=0;i<shell.ngp;i++)
  {
	*(L+i)=(double *)malloc(3*sizeof(double));
  }
  aa=(double *)malloc(3*sizeof(double));
  bb=(double *)malloc(3*sizeof(double));
  cc=(double *)malloc(3*sizeof(double));
  dd=(double *)malloc(3*sizeof(double));
  ee=(double *)malloc(3*sizeof(double));



  for(i=0;i<3;i++)
  {
	*(*(L+0)+i)=1.0/3.0;
	for(j=0;j<3;j++)
	{
	  if(i==j)
	  {
		*(*(L+(4+j))+i)=1.0;
	  }
	  else
	  {
		*(*(L+(4+j))+i)=0.0;
	  }

	  if((i+1)%3==j)
	  {
		*(*(L+(1+j))+i)=0.0;
	  }
	  else
	  {
		*(*(L+(1+j))+i)=0.5;
	  }
	}
  }
  /*area coordinate Li=(ai+bix+ciy)/(2*area)*/
  /*3 node triangular element with 7 integrated (Gauss) point*/
  /*
  *(*(L+0)+0)=1.0/3.0;*(*(L+0)+1)=1.0/3.0;*(*(L+0)+2)=1.0/3.0;
  *(*(L+1)+0)=1.0/2.0;*(*(L+1)+1)=1.0/2.0;*(*(L+1)+2)=0.0;
  *(*(L+2)+0)=0.0;*(*(L+2)+1)=1.0/2.0;*(*(L+2)+2)=1.0/2.0;
  *(*(L+3)+0)=1.0/2.0;*(*(L+3)+1)=0.0;*(*(L+3)+2)=1.0/2.0;
  *(*(L+4)+0)=1.0;*(*(L+4)+1)=0.0;*(*(L+4)+2)=0.0;
  *(*(L+5)+0)=0.0;*(*(Lz+5)+1)=1.0;*(*(L+5)+2)=0.0;
  *(*(L+6)+0)=0.0;*(*(L+6)+1)=0.0;*(*(L+6)+2)=1.0;
  */

  for(i=0;i<3;i++)/*NODE 0-0, 0-1, 0-2*/
  {
	for(j=0;j<2;j++)/*IN-PLANE 2D VECTOR*/
	{
	  *(*(exy+i)+j)=(shell.node[i]->d[0]-shell.node[0]->d[0])**(*(drccos+j)+0)
				   +(shell.node[i]->d[1]-shell.node[0]->d[1])**(*(drccos+j)+1)
				   +(shell.node[i]->d[2]-shell.node[0]->d[2])**(*(drccos+j)+2);
	  /*shell local coordinate {xi,yi(,zi=0)}*/
	}
  }
  area = shell.area;

  for(i=0;i<3;i++)
  {
	  j=i+1;
	  k=i+2;
	  if(j>=3)j-=3;
	  if(k>=3)k-=3;

	  *(a+i)=*(*(exy+j)+0)**(*(exy+k)+1)-*(*(exy+k)+0)**(*(exy+j)+1);/*xjyk-xkyj*/
	  *(b+i)=*(*(exy+j)+1)-*(*(exy+k)+1);/*yj-yk*/
	  *(c+i)=*(*(exy+k)+0)-*(*(exy+j)+0);/*xk-xj*/

	  /*CHECK : definition of A*/
	  *(Lx+i)=0.5**(b+i)/area;/*dLi/dx=bi/(2*area)*/
	  *(Ly+i)=0.5**(c+i)/area;/*dLi/dy=ci/(2*area)*/
	  /*area coordinate Li=(ai+bix+ciy)/(2*area)*/
  }

#if 1
  for(i=0;i<3;i++)
  {
	len    = *(b+i)**(b+i) + *(c+i)**(c+i);
	*(aa+i)= *(c+i)                                /len;
	*(bb+i)=-0.75**(b+i)**(c+i)                    /len;
	*(cc+i)=(0.25**(c+i)**(c+i)-0.5**(b+i)**(b+i)) /len;
	*(dd+i)=-*(b+i)                                /len;
	*(ee+i)=(0.25**(b+i)**(b+i)-0.5**(c+i)**(c+i)) /len;
  }
#endif


  Bp=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bp+i)=(double *)malloc(6*sizeof(double));
  }
  Bb=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bb+i)=(double *)malloc(9*sizeof(double));
  }
  Bt=(double **)malloc(1*sizeof(double *));
  for(i=0;i<1;i++)
  {
	*(Bt+i)=(double *)malloc(3*sizeof(double));
  }

  B=(double ***)malloc(shell.ngp*sizeof(double **));
  for(ii=0;ii<shell.ngp;ii++)
  {
	*(B+ii)=(double **)malloc(nstress*sizeof(double*));
	for(i=0;i<nstress;i++)
	{
	  *(*(B+ii)+i)=(double *)malloc(18*sizeof(double));
	  for(j=0;j<18;j++)
	  {
		*(*(*(B+ii)+i)+j) = 0.0;
      }
	}
  }

  /*SHAPE FUNCTION FOR IN-PLANE*/
  *(*(Bp+0)+0)=*(Lx+0);
  *(*(Bp+0)+1)=0.0;
  *(*(Bp+0)+2)=*(Lx+1);
  *(*(Bp+0)+3)=0.0;
  *(*(Bp+0)+4)=*(Lx+2);
  *(*(Bp+0)+5)=0.0;

  *(*(Bp+1)+0)=0.0;
  *(*(Bp+1)+1)=*(Ly+0);
  *(*(Bp+1)+2)=0.0;
  *(*(Bp+1)+3)=*(Ly+1);
  *(*(Bp+1)+4)=0.0;
  *(*(Bp+1)+5)=*(Ly+2);

  *(*(Bp+2)+0)=*(Ly+0);
  *(*(Bp+2)+1)=*(Lx+0);
  *(*(Bp+2)+2)=*(Ly+1);
  *(*(Bp+2)+3)=*(Lx+1);
  *(*(Bp+2)+4)=*(Ly+2);
  *(*(Bp+2)+5)=*(Lx+2);

  /*SHAPE FUNCTION FOR BENDING*/
  for(ii=0;ii<shell.ngp;ii++)
  {

#if 0
	/*Zienkiewicz‚ç‚Ì”ñ“K‡ŽOŠpŒ`—v‘fT-9N*/
	for(i=0;i<3;i++)
	{
	  j=(i+1)%3;
	  k=(i+2)%3;

	  *(*(Bb+0)+3*i+0)=-(
						1.0 * *(Lx+i)**(Lx+i) *   2.0*( *(*(L+ii)+j)+*(*(L+ii)+k) )
					   -1.0 * *(Lx+j)**(Lx+j) *   2.0*( *(*(L+ii)+i)              )
					   -1.0 * *(Lx+k)**(Lx+k) *   2.0*( *(*(L+ii)+i)              )
					   +2.0 * *(Lx+i)**(Lx+j) *   2.0*( *(*(L+ii)+i)-*(*(L+ii)+j) )
					   +0.0
					   +2.0 * *(Lx+k)**(Lx+i) *   2.0*( *(*(L+ii)+i)-*(*(L+ii)+k) )
					   );

	  *(*(Bb+0)+3*i+1)=-(
						1.0 * *(Lx+i)**(Lx+i) *   2.0*( *(b+j)**(*(L+ii)+k) - *(b+k)**(*(L+ii)+j) )
					   +0.0
					   +0.0
					   +2.0 * *(Lx+i)**(Lx+j) * (-2.0**(b+k)**(*(L+ii)+i) + 0.5*(*(b+j)-*(b+k))**(*(L+ii)+k) )
					   +2.0 * *(Lx+j)**(Lx+k) * ( 0.5*(*(b+j)-*(b+k))**(*(L+ii)+i) )
					   +2.0 * *(Lx+k)**(Lx+i) * ( 2.0**(b+j)**(*(L+ii)+i) + 0.5*(*(b+j)-*(b+k))**(*(L+ii)+j) )
					   );

	  *(*(Bb+0)+3*i+2)=-(
						1.0	* *(Lx+i)**(Lx+i) *   2.0*( *(c+j)**(*(L+ii)+k) - *(c+k)**(*(L+ii)+j) )
					   +0.0
					   +0.0
					   +2.0 * *(Lx+i)**(Lx+j) * (-2.0**(c+k)**(*(L+ii)+i) + 0.5*(*(c+j)-*(c+k))**(*(L+ii)+k) )
					   +2.0 * *(Lx+j)**(Lx+k) * ( 0.5*(*(c+j)-*(c+k))**(*(L+ii)+i) )
					   +2.0 * *(Lx+k)**(Lx+i) * ( 2.0**(c+j)**(*(L+ii)+i) + 0.5*(*(c+j)-*(c+k))**(*(L+ii)+j) )
					   );

	  *(*(Bb+1)+3*i+0)=-(
						1.0 * *(Ly+i)**(Ly+i) *   2.0*( *(*(L+ii)+j)+*(*(L+ii)+k) )
					   -1.0 * *(Ly+j)**(Ly+j) *   2.0*( *(*(L+ii)+i)              )
					   -1.0 * *(Ly+k)**(Ly+k) *   2.0*( *(*(L+ii)+i)              )
					   +2.0 * *(Ly+i)**(Ly+j) *   2.0*( *(*(L+ii)+i)-*(*(L+ii)+j) )
					   +0.0
					   +2.0 * *(Ly+k)**(Ly+i) *   2.0*( *(*(L+ii)+i)-*(*(L+ii)+k) )
					   );

	  *(*(Bb+1)+3*i+1)=-(
						1.0 * *(Ly+i)**(Ly+i) *   2.0*( *(b+j)**(*(L+ii)+k) - *(b+k)**(*(L+ii)+j) )
					   +0.0
					   +0.0
					   +2.0 * *(Ly+i)**(Ly+j) * (-2.0**(b+k)**(*(L+ii)+i) + 0.5*(*(b+j)-*(b+k))**(*(L+ii)+k) )
					   +2.0 * *(Ly+j)**(Ly+k) * ( 0.5*(*(b+j)-*(b+k))**(*(L+ii)+i) )
					   +2.0 * *(Ly+k)**(Ly+i) * ( 2.0**(b+j)**(*(L+ii)+i) + 0.5*(*(b+j)-*(b+k))**(*(L+ii)+j) )
					   );

	  *(*(Bb+1)+3*i+2)=-(
						1.0 * *(Ly+i)**(Ly+i) *   2.0*( *(c+j)**(*(L+ii)+k) -*(c+k)**(*(L+ii)+j) )
					   +0.0
					   +0.0
					   +2.0 * *(Ly+i)**(Ly+j) * (-2.0**(c+k)**(*(L+ii)+i) + 0.5*(*(c+j)-*(c+k))**(*(L+ii)+k) )
					   +2.0 * *(Ly+j)**(Ly+k) * ( 0.5*(*(c+j)-*(c+k))**(*(L+ii)+i) )
					   +2.0 * *(Ly+k)**(Ly+i) * ( 2.0**(c+j)**(*(L+ii)+i) + 0.5*(*(c+j)-*(c+k))**(*(L+ii)+j) )
					   );

	  *(*(Bb+2)+3*i+0)=-2.0*(
						   *(Lx+i)**(Ly+i)                     *  2.0*( *(*(L+ii)+j)+*(*(L+ii)+k) )
					   -   *(Lx+j)**(Ly+j)                     *  2.0*( *(*(L+ii)+i)              )
					   -   *(Lx+k)**(Ly+k)                     *  2.0*( *(*(L+ii)+i)              )
					   + ( *(Lx+i)**(Ly+j) + *(Lx+j)**(Ly+i) ) *  2.0*( *(*(L+ii)+i)-*(*(L+ii)+j) )
					   +   0.0
					   + ( *(Lx+k)**(Ly+i) + *(Lx+i)**(Ly+k) ) *  2.0*( *(*(L+ii)+i)-*(*(L+ii)+k) )
					   );

	  *(*(Bb+2)+3*i+1)=-2.0*(
						   *(Lx+i)**(Ly+i)                     *  2.0*( *(b+j)**(*(L+ii)+k) - *(b+k)**(*(L+ii)+j) )
					   +   0.0
					   +   0.0
					   + ( *(Lx+i)**(Ly+j) + *(Lx+j)**(Ly+i) ) *(-2.0**(b+k)**(*(L+ii)+i) + 0.5*(*(b+j)-*(b+k))**(*(L+ii)+k) )
					   + ( *(Lx+j)**(Ly+k) + *(Lx+k)**(Ly+j) ) *( 0.5*(*(b+j)-*(b+k))**(*(L+ii)+i) )
					   + ( *(Lx+k)**(Ly+i) + *(Lx+i)**(Ly+k) ) *( 2.0**(b+j)**(*(L+ii)+i) + 0.5*(*(b+j)-*(b+k))**(*(L+ii)+j) )
					   );

	  *(*(Bb+2)+3*i+2)=-2.0*(
						   *(Lx+i)**(Ly+i)                     *  2.0*( *(c+j)**(*(L+ii)+k) -*(c+k)**(*(L+ii)+j) )
					   +   0.0
					   +   0.0
					   + ( *(Lx+i)**(Ly+j) + *(Lx+j)**(Ly+i) ) *(-2.0**(c+k)**(*(L+ii)+i) + 0.5*(*(c+j)-*(c+k))**(*(L+ii)+k) )
					   + ( *(Lx+j)**(Ly+k) + *(Lx+k)**(Ly+j) ) *( 0.5*(*(c+j)-*(c+k))**(*(L+ii)+i) )
					   + ( *(Lx+k)**(Ly+i) + *(Lx+i)**(Ly+k) ) *( 2.0**(c+j)**(*(L+ii)+i) + 0.5*(*(c+j)-*(c+k))**(*(L+ii)+j) )
					   );
	}
#endif



#if 1
	/*Stricklin/Dhatt‚ç‚Ì—£ŽUKirchhoff‰¼’èŽOŠpŒ`(DKT:Discrete Kirchhoff Triangular)—v‘fT-9D*/
	for(i=0;i<3;i++)
	{
	  j=(i+1)%3;
	  k=(i+2)%3;

	  *(*(Bb+0)+3*i+0)=  6.0 * *(aa+k) * ( *(Lx+i) * *(*(L+ii)+j) + *(Lx+j) * *(*(L+ii)+i) )
					   - 6.0 * *(aa+j) * ( *(Lx+k) * *(*(L+ii)+i) + *(Lx+i) * *(*(L+ii)+k) );

	  *(*(Bb+0)+3*i+1)=  4.0 * *(bb+j) * ( *(Lx+k) * *(*(L+ii)+i) + *(Lx+i) * *(*(L+ii)+k) )
					   + 4.0 * *(bb+k) * ( *(Lx+i) * *(*(L+ii)+j) + *(Lx+j) * *(*(L+ii)+i) );

	  *(*(Bb+0)+3*i+2)=        *(Lx+i) * ( 4.0 * *(*(L+ii)+i) - 1.0 )
					   - 4.0 * *(cc+j) * ( *(Lx+k) * *(*(L+ii)+i) + *(Lx+i) * *(*(L+ii)+k) )
					   - 4.0 * *(cc+k) * ( *(Lx+i) * *(*(L+ii)+j) + *(Lx+j) * *(*(L+ii)+i) );

	  *(*(Bb+1)+3*i+0)=  6.0 * *(dd+k) * ( *(Ly+i) * *(*(L+ii)+j) + *(Ly+j) * *(*(L+ii)+i) )
					   - 6.0 * *(dd+j) * ( *(Ly+k) * *(*(L+ii)+i) + *(Ly+i) * *(*(L+ii)+k) );

	  *(*(Bb+1)+3*i+1)=        *(Ly+i) * ( 1.0 - 4.0 * *(*(L+ii)+i) )
					   + 4.0 * *(ee+j) * ( *(Ly+k) * *(*(L+ii)+i) + *(Ly+i) * *(*(L+ii)+k) )
					   + 4.0 * *(ee+k) * ( *(Ly+i) * *(*(L+ii)+j) + *(Ly+j) * *(*(L+ii)+i) );

	  *(*(Bb+1)+3*i+2)=- 4.0 * *(bb+j) * ( *(Ly+k) * *(*(L+ii)+i) + *(Ly+i) * *(*(L+ii)+k) )
					   - 4.0 * *(bb+k) * ( *(Ly+i) * *(*(L+ii)+j) + *(Ly+j) * *(*(L+ii)+i) );

	  *(*(Bb+2)+3*i+0)=  6.0 * *(aa+k) * ( *(Ly+i) * *(*(L+ii)+j) + *(Ly+j) * *(*(L+ii)+i) )
					   - 6.0 * *(aa+j) * ( *(Ly+k) * *(*(L+ii)+i) + *(Ly+i) * *(*(L+ii)+k) )
					   + 6.0 * *(dd+k) * ( *(Lx+i) * *(*(L+ii)+j) + *(Lx+j) * *(*(L+ii)+i) )
					   - 6.0 * *(dd+j) * ( *(Lx+k) * *(*(L+ii)+i) + *(Lx+i) * *(*(L+ii)+k) );

	  *(*(Bb+2)+3*i+1)=  4.0 * *(bb+j) * ( *(Ly+k) * *(*(L+ii)+i) + *(Ly+i) * *(*(L+ii)+k) )
					   + 4.0 * *(bb+k) * ( *(Ly+i) * *(*(L+ii)+j) + *(Ly+j) * *(*(L+ii)+i) )
					   +       *(Lx+i) * ( 1.0 - 4.0 * *(*(L+ii)+i) )
					   + 4.0 * *(ee+j) * ( *(Lx+k) * *(*(L+ii)+i) + *(Lx+i) * *(*(L+ii)+k) )
					   + 4.0 * *(ee+k) * ( *(Lx+i) * *(*(L+ii)+j) + *(Lx+j) * *(*(L+ii)+i) );

	  *(*(Bb+2)+3*i+2)=        *(Ly+i) * ( 4.0 * *(*(L+ii)+i) - 1.0 )
					   - 4.0 * *(cc+j) * ( *(Ly+k) * *(*(L+ii)+i) + *(Ly+i) * *(*(L+ii)+k) )
					   - 4.0 * *(cc+k) * ( *(Ly+i) * *(*(L+ii)+j) + *(Ly+j) * *(*(L+ii)+i) )
					   - 4.0 * *(bb+j) * ( *(Lx+k) * *(*(L+ii)+i) + *(Lx+i) * *(*(L+ii)+k) )
					   - 4.0 * *(bb+k) * ( *(Lx+i) * *(*(L+ii)+j) + *(Lx+j) * *(*(L+ii)+i) );
	}
#endif

	  *(*(Bt+0)+0)=*(*(L+ii)+0)-*(*(L+0)+0);
	  *(*(Bt+0)+1)=*(*(L+ii)+1)-*(*(L+0)+1);
	  *(*(Bt+0)+2)=*(*(L+ii)+2)-*(*(L+0)+2);


	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(*(B+ii)+i)+6*j+0)=*(*(Bp+i)+2*j+0);
		*(*(*(B+ii)+i)+6*j+1)=*(*(Bp+i)+2*j+1);

		*(*(*(B+ii)+3+i)+6*j+2)=*(*(Bb+i)+3*j+0);
		*(*(*(B+ii)+3+i)+6*j+3)=*(*(Bb+i)+3*j+1);
		*(*(*(B+ii)+3+i)+6*j+4)=*(*(Bb+i)+3*j+2);
	  }
	}
	for(i=0;i<1;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(*(B+ii)+6+i)+6*j+5)=*(*(Bt+i)+j);
	  }
	}
  }


  free(a);
  free(b);
  free(c);
  free(Lx);
  free(Ly);
  freematrix(exy,3);
  freematrix(L,shell.ngp);
  free(aa);
  free(bb);
  free(cc);
  free(dd);
  free(ee);
  freematrix(Bp,3);
  freematrix(Bb,3);
  freematrix(Bt,1);

  return B;
}


double **assemshellemtx(struct oshell shell)
/*ASSEMBLAGE ELASTIC MATRIX.*/
{
  int i,j;
  int ngp = shell.ngp;
  int nstress = shell.nstress;
  int ii;/*LOOP FOR INTEGRATION POINTS*/
  double** Ke;
  double area;
  double* w;
  double** drccos;
  double*** C, ***B;
  double** Dp,** Bp,** Bpt,** DpBp,** Kp;
  double** Db,** Bb,** Bbt,** DbBb,** Kb;
  double** Dt,** Bt,** Btt,** DtBt,** Kt;

  double prate = 1.0/*shell.prate*/;
  double brate = 1.0/*shell.brate*/;

  Ke=(double **)malloc(18*sizeof(double *));
  for(i=0;i<18;i++)
  {
	*(Ke+i)=(double *)malloc(18*sizeof(double));
	for(j=0;j<18;j++)
	{
	  *(*(Ke+i)+j)=0.0;                                              /*INITIAL.*/
	}
  }

  area = shell.area;
  w=(double *)malloc(ngp*sizeof(double));
  for(ii=0;ii<ngp;ii++)
  {
	*(w+ii)=area*shell.w[ii];
  }



  Dp=(double **)malloc(3*sizeof(double *));
  Db=(double **)malloc(3*sizeof(double *));
  Dt=(double **)malloc(1*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Dp+i)=(double *)malloc(3*sizeof(double));
  }
  for(i=0;i<3;i++)
  {
	*(Db+i)=(double *)malloc(3*sizeof(double));
  }
  for(i=0;i<1;i++)
  {
	*(Dt+i)=(double *)malloc(1*sizeof(double));
  }

  Bp=(double **)malloc(3*sizeof(double *));
  Bb=(double **)malloc(3*sizeof(double *));
  Bt=(double **)malloc(1*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bp+i)=(double *)malloc(6*sizeof(double));
  }
  for(i=0;i<3;i++)
  {
	*(Bb+i)=(double *)malloc(9*sizeof(double));
  }
  for(i=0;i<1;i++)
  {
	*(Bt+i)=(double *)malloc(3*sizeof(double));
  }

  drccos = shelldrccos(shell);
  C = elasticCshell(shell);
  B = assemshellshape(shell, drccos);






	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Dp+i)+j) = *(*(*(C+0)+i)+j);
		*(*(Db+i)+j) = *(*(*(C+0)+3+i)+3+j);
	  }
	}
	*(*(Dt+0)+0) = *(*(*(C+0)+6)+6);

	//IN-PLANE
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bp+i)+2*j+0) = *(*(*(B+0)+i)+6*j+0);
		*(*(Bp+i)+2*j+1) = *(*(*(B+0)+i)+6*j+1);
	  }
	}

	Bpt=matrixtransposeIII(Bp,3,6);
	DpBp=matrixmatrixIII(Dp,Bp,3,3,6);
	Kp=matrixmatrixIII(Bpt,DpBp,6,3,6);


	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Ke+6*i+0)+6*j+0)+=*(*(Kp+2*i+0)+2*j+0)*area*prate;
		*(*(Ke+6*i+0)+6*j+1)+=*(*(Kp+2*i+0)+2*j+1)*area*prate;
		*(*(Ke+6*i+1)+6*j+0)+=*(*(Kp+2*i+1)+2*j+0)*area*prate;
		*(*(Ke+6*i+1)+6*j+1)+=*(*(Kp+2*i+1)+2*j+1)*area*prate;
	  }
	}

	freematrix(Bpt,6);
	freematrix(DpBp,3);
	freematrix(Kp,6);


  for(ii=0;ii<ngp;ii++)
  {
	//BENDING
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bb+i)+3*j+0) = *(*(*(B+ii)+3+i)+6*j+2);
		*(*(Bb+i)+3*j+1) = *(*(*(B+ii)+3+i)+6*j+3);
		*(*(Bb+i)+3*j+2) = *(*(*(B+ii)+3+i)+6*j+4);
	  }
	}

	Bbt=matrixtransposeIII(Bb,3,9);
	DbBb=matrixmatrixIII(Db,Bb,3,3,9);
	Kb=matrixmatrixIII(Bbt,DbBb,9,3,9);

	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Ke+6*i+2)+6*j+2) += *(*(Kb+3*i+0)+3*j+0)**(w+ii)*brate;
		*(*(Ke+6*i+2)+6*j+3) += *(*(Kb+3*i+0)+3*j+1)**(w+ii)*brate;
		*(*(Ke+6*i+2)+6*j+4) += *(*(Kb+3*i+0)+3*j+2)**(w+ii)*brate;
		*(*(Ke+6*i+3)+6*j+2) += *(*(Kb+3*i+1)+3*j+0)**(w+ii)*brate;
		*(*(Ke+6*i+3)+6*j+3) += *(*(Kb+3*i+1)+3*j+1)**(w+ii)*brate;
		*(*(Ke+6*i+3)+6*j+4) += *(*(Kb+3*i+1)+3*j+2)**(w+ii)*brate;
		*(*(Ke+6*i+4)+6*j+2) += *(*(Kb+3*i+2)+3*j+0)**(w+ii)*brate;
		*(*(Ke+6*i+4)+6*j+3) += *(*(Kb+3*i+2)+3*j+1)**(w+ii)*brate;
		*(*(Ke+6*i+4)+6*j+4) += *(*(Kb+3*i+2)+3*j+2)**(w+ii)*brate;
	  }
	}

	freematrix(Bbt,9);
	freematrix(DbBb,3);
	freematrix(Kb,9);

	//TORSION
	for(i=0;i<1;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bt+i)+j) = *(*(*(B+ii)+6+i)+6*j+5);
	  }
	}

	Btt=matrixtransposeIII(Bt,1,3);
	DtBt=matrixmatrixIII(Dt,Bt,1,1,3);
	Kt=matrixmatrixIII(Btt,DtBt,3,1,3);

	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Ke+6*i+5)+6*j+5) += *(*(Kt+i)+j)**(w+ii);
	  }
	}

	freematrix(Btt,3);
	freematrix(DtBt,1);
	freematrix(Kt,3);

  }

  /*DRILL*/

  /*
  double E=shell.sect->E;
  double t=shell.sect->area;

  for(ii=0;ii<1;ii++)
  {
	*(*(Ke+ 5)+ 5)= 0.030*E*t*area;
	*(*(Ke+ 5)+11)=-0.015*E*t*area;
	*(*(Ke+ 5)+17)=*(*(Ke+5)+11);
	*(*(Ke+11)+ 5)=*(*(Ke+5)+11);
	*(*(Ke+11)+11)=*(*(Ke+5)+ 5);
	*(*(Ke+11)+17)=*(*(Ke+5)+11);
	*(*(Ke+17)+ 5)=*(*(Ke+5)+11);
	*(*(Ke+17)+11)=*(*(Ke+5)+11);
	*(*(Ke+17)+17)=*(*(Ke+5)+ 5);
  }
  */


  freematrix(drccos,3);
  for(ii=0;ii<ngp;ii++)freematrix(*(C+ii), nstress);
  free(C);
  for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
  free(B);

  free(w);
  freematrix(Dp,3);
  freematrix(Db,3);
  freematrix(Dt,1);
  freematrix(Bp,3);
  freematrix(Bb,3);
  freematrix(Bt,1);

  return Ke;
}


double **assemshellpmtx(struct oshell shell,double*** C,double *** B)
/*ASSEMBLAGE ELASTO-PLASTIC MATRIX.*/
{
  int i,j;
  int ii;/*LOOP FOR INTEGRATION POINTS*/
  int ngp = shell.ngp;
  int nstress = shell.nstress;
  double** K;
  double* w;
  double area;
  double** Dp,** Bp,** Bpt,** DpBp,** Kp;
  double** Db,** Bb,** Bbt,** DbBb,** Kb;
  double** Dt,** Bt,** Btt,** DtBt,** Kt;
  double** Bat,**DaBa,**Ka;

  double prate = shell.prate;
  double brate = shell.brate;

  K=(double **)malloc(18*sizeof(double *));
  for(i=0;i<18;i++)
  {
	*(K+i)=(double *)malloc(18*sizeof(double));
	for(j=0;j<18;j++)
	{
	  *(*(K+i)+j)=0.0;                                              /*INITIAL.*/
	}
  }

  area = shell.area;
  w=(double *)malloc(ngp*sizeof(double));
  for(ii=0;ii<ngp;ii++)
  {
	*(w+ii)=area*shell.w[ii];
  }


  Dp=(double **)malloc(3*sizeof(double *));
  Db=(double **)malloc(3*sizeof(double *));
  Dt=(double **)malloc(1*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Dp+i)=(double *)malloc(3*sizeof(double));
  }
  for(i=0;i<3;i++)
  {
	*(Db+i)=(double *)malloc(3*sizeof(double));
  }
  for(i=0;i<1;i++)
  {
	*(Dt+i)=(double *)malloc(1*sizeof(double));
  }

  Bp=(double **)malloc(3*sizeof(double *));
  Bb=(double **)malloc(3*sizeof(double *));
  Bt=(double **)malloc(1*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bp+i)=(double *)malloc(6*sizeof(double));
  }
  for(i=0;i<3;i++)
  {
	*(Bb+i)=(double *)malloc(9*sizeof(double));
  }
  for(i=0;i<1;i++)
  {
	*(Bt+i)=(double *)malloc(3*sizeof(double));
  }



  for(ii=0;ii<ngp;ii++)
  {
	/*
	//STRESS-STRAIN MATRIX
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Dp+i)+j) = *(*(*(C+ii)+i)+j);
		*(*(Db+i)+j) = *(*(*(C+ii)+3+i)+3+j);
	  }
	}
	*(*(Dt+0)+0) = *(*(*(C+ii)+6)+6);

	//IN-PLANE
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bp+i)+2*j+0) = *(*(*(B+ii)+i)+6*j+0);
		*(*(Bp+i)+2*j+1) = *(*(*(B+ii)+i)+6*j+1);
	  }
	}

	Bpt=matrixtransposeIII(Bp,3,6);
	DpBp=matrixmatrixIII(Dp,Bp,3,3,6);
	Kp=matrixmatrixIII(Bpt,DpBp,6,3,6);


	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(K+6*i+0)+6*j+0)+=*(*(Kp+2*i+0)+2*j+0)**(w+ii)*prate;
		*(*(K+6*i+0)+6*j+1)+=*(*(Kp+2*i+0)+2*j+1)**(w+ii)*prate;
		*(*(K+6*i+1)+6*j+0)+=*(*(Kp+2*i+1)+2*j+0)**(w+ii)*prate;
		*(*(K+6*i+1)+6*j+1)+=*(*(Kp+2*i+1)+2*j+1)**(w+ii)*prate;
	  }
	}

	freematrix(Bpt,6);
	freematrix(DpBp,3);
	freematrix(Kp,6);

	//BENDING
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bb+i)+3*j+0) = *(*(*(B+ii)+3+i)+6*j+2);
		*(*(Bb+i)+3*j+1) = *(*(*(B+ii)+3+i)+6*j+3);
		*(*(Bb+i)+3*j+2) = *(*(*(B+ii)+3+i)+6*j+4);
	  }
	}

	Bbt=matrixtransposeIII(Bb,3,9);
	DbBb=matrixmatrixIII(Db,Bb,3,3,9);
	Kb=matrixmatrixIII(Bbt,DbBb,9,3,9);

	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(K+6*i+2)+6*j+2) += *(*(Kb+3*i+0)+3*j+0)**(w+ii)*brate;
		*(*(K+6*i+2)+6*j+3) += *(*(Kb+3*i+0)+3*j+1)**(w+ii)*brate;
		*(*(K+6*i+2)+6*j+4) += *(*(Kb+3*i+0)+3*j+2)**(w+ii)*brate;
		*(*(K+6*i+3)+6*j+2) += *(*(Kb+3*i+1)+3*j+0)**(w+ii)*brate;
		*(*(K+6*i+3)+6*j+3) += *(*(Kb+3*i+1)+3*j+1)**(w+ii)*brate;
		*(*(K+6*i+3)+6*j+4) += *(*(Kb+3*i+1)+3*j+2)**(w+ii)*brate;
		*(*(K+6*i+4)+6*j+2) += *(*(Kb+3*i+2)+3*j+0)**(w+ii)*brate;
		*(*(K+6*i+4)+6*j+3) += *(*(Kb+3*i+2)+3*j+1)**(w+ii)*brate;
		*(*(K+6*i+4)+6*j+4) += *(*(Kb+3*i+2)+3*j+2)**(w+ii)*brate;
	  }
	}

	freematrix(Bbt,9);
	freematrix(DbBb,3);
	freematrix(Kb,9);

	//TORSION
	for(i=0;i<1;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bt+i)+j) = *(*(*(B+ii)+6+i)+6*j+5);
	  }
	}

	Btt=matrixtransposeIII(Bt,1,3);
	DtBt=matrixmatrixIII(Dt,Bt,1,1,3);
	Kt=matrixmatrixIII(Btt,DtBt,3,1,3);

	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(K+6*i+5)+6*j+5) += *(*(Kt+i)+j)**(w+ii);
	  }
	}

	freematrix(Btt,3);
	freematrix(DtBt,1);
	freematrix(Kt,3);

	*/

	Bat=matrixtransposeIII(*(B+ii),7,18);
	DaBa=matrixmatrixIII(*(C+ii),*(B+ii),7,7,18);
	Ka=matrixmatrixIII(Bat,DaBa,18,7,18);

	for(i=0;i<18;i++)
	{
	  for(j=0;j<18;j++)
	  {
		*(*(K+i)+j) += *(*(Ka+i)+j)**(w+ii);
	  }
	}
	freematrix(Bat,18);
	freematrix(DaBa,7);
	freematrix(Ka,18);




  }


  free(w);
  freematrix(Dp,3);
  freematrix(Db,3);
  freematrix(Dt,1);
  freematrix(Bp,3);
  freematrix(Bb,3);
  freematrix(Bt,1);

  return K;
}

double **assemshellmmtx(struct oshell shell,double **drccos)
/*ASSEMBLAGE ELASTIC MATRIX.*/
{
  char string[100];
  int i,j,k,ii;
  double **m,**exy,**N,**Nt,**NtN,**L;
  double t,hiju;
  double *b,*c;
  double *a;
  double det;
  double Liij,Lijj,Liik,Likk,Lijk;
  double alpha=0.0;
  //double alpha=0.15;

  t=shell.sect->area;
  hiju=shell.sect->hiju[0];

  m=(double **)malloc(18*sizeof(double *));
  for(i=0;i<18;i++)
  {
	*(m+i)=(double *)malloc(18*sizeof(double));
	for(j=0;j<18;j++)
	{
	  *(*(m+i)+j)=0.0;                                              /*INITIAL.*/
	}
  }

  exy=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(exy+i)=(double *)malloc(2*sizeof(double));
	for(j=0;j<2;j++)
	{
	  *(*(exy+i)+j)=(shell.node[i]->d[0]-shell.node[0]->d[0])**(*(drccos+j)+0)
				   +(shell.node[i]->d[1]-shell.node[0]->d[1])**(*(drccos+j)+1)
				   +(shell.node[i]->d[2]-shell.node[0]->d[2])**(*(drccos+j)+2);
	}
  }
  /*shell local coordinate {xi,yi(,zi=0)}*/
  det=0.5*(*(*(exy+1)+0)**(*(exy+2)+1)-*(*(exy+1)+1)**(*(exy+2)+0));

#if 1
  /*WEIGHT OF EACH NODE*/
  a=(double *)malloc(7*sizeof(double));
  *(a+0)=27.0/60.0;
  *(a+1)=8.0/60.0;
  *(a+2)=8.0/60.0;
  *(a+3)=8.0/60.0;
  *(a+4)=3.0/60.0;
  *(a+5)=3.0/60.0;
  *(a+6)=3.0/60.0;

  /*TRIANGLE COORDINATION*/
  L=(double **)malloc(7*sizeof(double *));
  for(i=0;i<7;i++)
  {
	*(L+i)=(double *)malloc(3*sizeof(double));
  }
  for(i=0;i<3;i++)
  {
	*(*(L+0)+i)=1.0/3.0;
	for(j=0;j<3;j++)
	{
	  if(i==j)
	  {
		*(*(L+(4+j))+i)=1.0;
	  }
	  else
	  {
		*(*(L+(4+j))+i)=0.0;
	  }

	  if((i+1)%3==j)
	  {
		*(*(L+(1+j))+i)=0.0;
	  }
	  else
	  {
		*(*(L+(1+j))+i)=0.5;
	  }
	}
  }
  /*area coordinate Li=(ai+bix+ciy)/(2*det)*/

  b=(double *)malloc(3*sizeof(double));
  c=(double *)malloc(3*sizeof(double));
  for(i=0;i<3;i++)
  {
	  j=i+1;
	  k=i+2;
	  if(j>=3)j-=3;
	  if(k>=3)k-=3;
	  *(b+i)=*(*(exy+j)+1)-*(*(exy+k)+1);/*yj-yk*/
	  *(c+i)=*(*(exy+k)+0)-*(*(exy+j)+0);/*xk-xj*/
  }



  for(ii=0;ii<7;ii++)
  {
	N=(double **)malloc(3*sizeof(double *));
	for(i=0;i<3;i++)
	{
	  *(N+i)=(double *)malloc(18*sizeof(double));
	  for(j=0;j<18;j++)
	  {
		*(*(N+i)+j)=0.0;                                              /*INITIAL.*/
	  }
	}
	for(i=0;i<3;i++)
	{
		j=i+1;
		k=i+2;
		if(j>=3)j-=3;
		if(k>=3)k-=3;
		Liij=*(*(L+ii)+i)**(*(L+ii)+i)**(*(L+ii)+j);
		Lijj=*(*(L+ii)+i)**(*(L+ii)+j)**(*(L+ii)+j);
		Liik=*(*(L+ii)+i)**(*(L+ii)+i)**(*(L+ii)+k);
		Likk=*(*(L+ii)+i)**(*(L+ii)+k)**(*(L+ii)+k);
		Lijk=*(*(L+ii)+i)**(*(L+ii)+j)**(*(L+ii)+k);

		*(*(N+0)+6*i+0) = *(*(L+ii)+i);
		*(*(N+1)+6*i+1) = *(*(L+ii)+i);
		*(*(N+2)+6*i+2) = *(*(L+ii)+i) + (1.0-alpha)*(Liij - Lijj + Liik - Likk);
		*(*(N+2)+6*i+3) = (1.0-alpha)*(*(b+j) * (Liik+0.5*Lijk) - *(b+k) * (Liij+0.5*Lijk));
		*(*(N+2)+6*i+4) = (1.0-alpha)*(*(c+j) * (Liik+0.5*Lijk) - *(c+k) * (Liij+0.5*Lijk));
	}
	Nt=matrixtransposeIII(N,3,18);
	NtN=matrixmatrixIII(Nt,N,18,3,18);
	for(i=0;i<18;i++)
	{
	  for(j=0;j<18;j++)
	  {
		*(*(m+i)+j) += *(*(NtN+i)+j)**(a+ii)*det*t*hiju;
	  }
	}
	freematrix(N,3);
	freematrix(Nt,18);
	freematrix(NtN,18);
  }

  free(a);
  freematrix(L,7);
  free(b);
  free(c);
#endif
#if 0
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  for(k=0;k<3;k++)
	  {
		if(i==j)
		{
		  *(*(m+6*i+k+0)+6*j+k+0) = hiju*det*t/6.0;
		  *(*(m+6*i+k+3)+6*j+k+3) = hiju*det*pow(t,3)/18.0;
		}
		else
		{
		  *(*(m+6*i+k+0)+6*j+k+0) = hiju*det*t/12.0;
		  *(*(m+6*i+k+3)+6*j+k+3) = hiju*det*pow(t,3)/36.0;
		}
	  }
	}
  }
  /*MASS MATRIX BY ZHONG AND ALMEDIA.*/
#endif

  freematrix(exy,3);

  return m;
}/*assemshellmmtx*/


double *assemshellpvct(struct oshell shell,double **drccos)
{/*nodal force equivalent to surface force*/
  int i,j,k,ii;
  double *q,*p,**exy,**Nt,**L;
  double t;
  double *b,*c;
  double *a;
  double det;
  double Liij,Lijj,Liik,Likk,Lijk;
  double *perpl,*lload;

  a=(double *)malloc(7*sizeof(double));
  *(a+0)=27.0/60.0;
  *(a+1)=8.0/60.0;
  *(a+2)=8.0/60.0;
  *(a+3)=8.0/60.0;
  *(a+4)=3.0/60.0;
  *(a+5)=3.0/60.0;
  *(a+6)=3.0/60.0;

  L=(double **)malloc(7*sizeof(double *));
  for(i=0;i<7;i++)
  {
	*(L+i)=(double *)malloc(3*sizeof(double));
  }
  for(i=0;i<3;i++)
  {
	*(*(L+0)+i)=1.0/3.0;
	for(j=0;j<3;j++)
	{
	  if(i==j)
	  {
		*(*(L+(4+j))+i)=1.0;
	  }
	  else
	  {
		*(*(L+(4+j))+i)=0.0;
	  }

	  if((i+1)%3==j)
	  {
		*(*(L+(1+j))+i)=0.0;
	  }
	  else
	  {
		*(*(L+(1+j))+i)=0.5;
	  }
	}
  }


  b=(double *)malloc(3*sizeof(double));
  c=(double *)malloc(3*sizeof(double));
  q=(double *)malloc(3*sizeof(double));
  exy=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(exy+i)=(double *)malloc(2*sizeof(double));
	for(j=0;j<2;j++)
	{
	  *(*(exy+i)+j)=(shell.node[i]->d[0]-shell.node[0]->d[0])**(*(drccos+j)+0)
				   +(shell.node[i]->d[1]-shell.node[0]->d[1])**(*(drccos+j)+1)
				   +(shell.node[i]->d[2]-shell.node[0]->d[2])**(*(drccos+j)+2);
	}
  }
  det=0.5*(*(*(exy+1)+0)**(*(exy+2)+1)-*(*(exy+1)+1)**(*(exy+2)+0));

  for(i=0;i<3;i++)
  {
	  j=i+1;
	  k=i+2;
	  if(j>=3)j-=3;
	  if(k>=3)k-=3;
	  *(b+i)=*(*(exy+j)+1)-*(*(exy+k)+1);
	  *(c+i)=*(*(exy+k)+0)-*(*(exy+j)+0);
  }

  t=shell.sect->area;
  perpl=shell.sect->perpl;
  lload=shell.sect->lload;
  for(i=0;i<3;i++)
  {
	*(q+i)=*(perpl+i);
	for(j=0;j<3;j++)
	{
	  *(q+i)+=*(*(drccos+i)+j)**(lload+j);
	}
  }


  p=(double *)malloc(18*sizeof(double));
  for(i=0;i<18;i++)
  {
	*(p+i)=0.0;
  }

  for(ii=0;ii<7;ii++)
  {
	Nt=(double **)malloc(18*sizeof(double *));
	for(i=0;i<18;i++)
	{
	  *(Nt+i)=(double *)malloc(3*sizeof(double));
	  for(j=0;j<3;j++)
	  {
		*(*(Nt+i)+j)=0.0;                                              /*INITIAL.*/
	  }
	}
	for(i=0;i<3;i++)
	{
		j=i+1;
		k=i+2;
		if(j>=3)j-=3;
		if(k>=3)k-=3;
		Liij=*(*(L+ii)+i)**(*(L+ii)+i)**(*(L+ii)+j);
		Lijj=*(*(L+ii)+i)**(*(L+ii)+j)**(*(L+ii)+j);
		Liik=*(*(L+ii)+i)**(*(L+ii)+i)**(*(L+ii)+k);
		Likk=*(*(L+ii)+i)**(*(L+ii)+k)**(*(L+ii)+k);
		Lijk=*(*(L+ii)+i)**(*(L+ii)+j)**(*(L+ii)+k);

		*(*(Nt+6*i+0)+0) = *(*(L+ii)+i);
		*(*(Nt+6*i+1)+1) = *(*(L+ii)+i);
		*(*(Nt+6*i+2)+2) = *(*(L+ii)+i) + Liij - Lijj + Liik - Likk;
		*(*(Nt+6*i+3)+2) = *(b+j) * (Liik+0.5*Lijk) - *(b+k) * (Liij+0.5*Lijk);
		*(*(Nt+6*i+4)+2) = *(c+j) * (Liik+0.5*Lijk) - *(c+k) * (Liij+0.5*Lijk);
	}

	for(i=0;i<18;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(p+i) += *(*(Nt+i)+j)**(q+j)**(a+ii)*det;
	  }
	}
	freematrix(Nt,18);
  }
  freematrix(L,7);
  freematrix(exy,3);
  free(a);
  free(b);
  free(c);
  free(q);
  return p;
}/*assemshellpvct*/


double shellvolume(struct oshell shell, double** drccos)
{
	double volume;
	struct onode node = *(shell.node[0]);

	volume = (*(*(drccos + 2) + 0) * node.d[GX]
		   +  *(*(drccos + 2) + 1) * node.d[GY]
		   +  *(*(drccos + 2) + 2) * node.d[GZ]) * shell.area / 3.0;
	return volume;
}/*shellvolume*/

void assemshellvolume(struct oshell* shells, int nshell, double* ddisp, double* volume)
{
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	double* gforminit, * eforminit;                      /*DEFORMED COORDINATION OF ELEMENT*/
	double* gform, * eform;                      /*INITIAL COORDINATION OF ELEMENT*/
	double** drccos,** drccosinit;                           /*MATRIX*/
	double area;

	for (i = 1; i <= nshell; i++)
	{
		inputshell(shells, NULL, i - 1, &shell);
		nnod = shell.nnod;


		/*INITIAL CONFIGURATION*/
		/*
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(iform, shell.node[ii]);
		}
		drccosinit = shelldrccos(shell, &area);
		gforminit = extractshelldisplacement(shell, iform);
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);
		*/

		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell);
		//gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
		//eform = extractlocalcoord(gform,drccos,nnod); 			       	    /*{Xe+Ue}*/

		*volume += shellvolume(shell, drccos);                   		/*VOLUME*/

		freematrix(drccos, 3);
		/*
		freematrix(drccosinit, 3);

		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		*/

	}
	return;
}


void inputshell(struct oshell *shells,
				struct memoryshell *mshell,int offset,
				struct oshell *shell)
{
  int i,j;
  shell->nnod=(shells+offset)->nnod;
  shell->ngp=(shells+offset)->ngp;
  shell->nstress=(shells+offset)->nstress;
  shell->loff=offset;
  shell->code=(shells+offset)->code;                  /*ELEMENT CODE.*/
  shell->sect=(shells+offset)->sect;               /*SECTION POINTER.*/
  shell->prate=(shells+offset)->prate;
  shell->brate=(shells+offset)->brate;
  shell->area=(shells+offset)->area;

  for(i=0;i<shell->nnod;i++)
  {
	shell->node[i]=(shells+offset)->node[i];

	if(mshell!=NULL)
	{
	  for(j=0;j<6;j++)                                      /*STRESS.*/
	  {
		shell->stress[i][j]=(mshell+offset)->stress[i][j];
	  }
	}
	else
	{
	  for(j=0;j<6;j++)                                      /*STRESS.*/
	  {
		shell->stress[i][j]=(shells+offset)->stress[i][j];
	  }
	}
  }
  for(i=0;i<shell->ngp;i++)
  {
	shell->w[i]=(shells+offset)->w[i];

	if(mshell!=NULL)
	{
	  for(j=0;j<shell->nstress;j++)                                      /*STRESS.*/
	  {
		(shell->gp[i]).estrain[j]=((mshell+offset)->gp[i]).estrain[j];
		(shell->gp[i]).pstrain[j]=((mshell+offset)->gp[i]).pstrain[j];
		(shell->gp[i]).stress[j]=((mshell+offset)->gp[i]).stress[j];
		(shell->gp[i]).backstress[j]=((mshell+offset)->gp[i]).backstress[j];
	  }
	  (shell->gp[i]).qn=((mshell+offset)->gp[i]).qn;
	  (shell->gp[i]).qm=((mshell+offset)->gp[i]).qm;
	  (shell->gp[i]).qnm=((mshell+offset)->gp[i]).qnm;
	  (shell->gp[i]).qilyushin=((mshell+offset)->gp[i]).qilyushin;
	  (shell->gp[i]).yinit=((mshell+offset)->gp[i]).yinit;
	  (shell->gp[i]).y=((mshell+offset)->gp[i]).y;

	  (shell->gp[i]).alpha=((mshell+offset)->gp[i]).alpha;
	  (shell->gp[i]).Ee=((mshell+offset)->gp[i]).Ee;
	  (shell->gp[i]).Ep=((mshell+offset)->gp[i]).Ep;
	}
	else
	{
	  for(j=0;j<shell->nstress;j++)                                      /*STRESS.*/
	  {
		(shell->gp[i]).estrain[j]=((shells+offset)->gp[i]).estrain[j];
		(shell->gp[i]).pstrain[j]=((shells+offset)->gp[i]).pstrain[j];
		(shell->gp[i]).stress[j]=((shells+offset)->gp[i]).stress[j];
		(shell->gp[i]).backstress[j]=((shells+offset)->gp[i]).backstress[j];
	  }
	  (shell->gp[i]).qn=((shells+offset)->gp[i]).qn;
	  (shell->gp[i]).qm=((shells+offset)->gp[i]).qm;
	  (shell->gp[i]).qnm=((shells+offset)->gp[i]).qnm;
	  (shell->gp[i]).qilyushin=((shells+offset)->gp[i]).qilyushin;
	  (shell->gp[i]).yinit=((shells+offset)->gp[i]).yinit;
	  (shell->gp[i]).y=((shells+offset)->gp[i]).y;

	  (shell->gp[i]).alpha=((shells+offset)->gp[i]).alpha;
	  (shell->gp[i]).Ee=((shells+offset)->gp[i]).Ee;
	  (shell->gp[i]).Ep=((shells+offset)->gp[i]).Ep;
	}
  }
  return;
}/*inputshell*/

void outputshell(struct oshell *shells,
				 int offset,
				 struct oshell *shell)
{
  int i,j;

  for(i=0;i<shell->nnod;i++)
  {
	(shells+offset)->node[i]=shell->node[i];
	for(j=0;j<6;j++)                                      /*STRESS.*/
	{
	  (shells+offset)->stress[i][j]=shell->stress[i][j];
	}
  }
  for(i=0;i<shell->ngp;i++)
  {
	for(j=0;j<shell->nstress;j++)                                      /*STRESS.*/
	{
	  ((shells+offset)->gp[i]).estrain[j]=(shell->gp[i]).estrain[j];
	  ((shells+offset)->gp[i]).pstrain[j]=(shell->gp[i]).pstrain[j];
	  ((shells+offset)->gp[i]).stress[j]=(shell->gp[i]).stress[j];
	  ((shells+offset)->gp[i]).backstress[j]=(shell->gp[i]).backstress[j];
	}
	((shells+offset)->gp[i]).qn=(shell->gp[i]).qn;
	((shells+offset)->gp[i]).qm=(shell->gp[i]).qm;
	((shells+offset)->gp[i]).qnm=(shell->gp[i]).qnm;
	((shells+offset)->gp[i]).qilyushin=(shell->gp[i]).qilyushin;
	((shells+offset)->gp[i]).yinit=(shell->gp[i]).yinit;
	((shells+offset)->gp[i]).y=(shell->gp[i]).y;

	((shells+offset)->gp[i]).alpha=(shell->gp[i]).alpha;
	((shells+offset)->gp[i]).Ee=(shell->gp[i]).Ee;
	((shells+offset)->gp[i]).Ep=(shell->gp[i]).Ep;
  }
  return;
}/*inputshell*/

void outputmemoryshell(struct oshell *shells,
					   struct memoryshell *mshell,int offset)
{
  int i,j;

  for(i=0;i<(shells+offset)->nnod;i++)
  {
	for(j=0;j<6;j++)                                      /*STRESS.*/
	{
	  (mshell+offset)->stress[i][j]=(shells+offset)->stress[i][j];
	}
  }
  for(i=0;i<(shells+offset)->ngp;i++)
  {
	for(j=0;j<(shells+offset)->nstress;j++)                                      /*STRESS.*/
	{
	  ((mshell+offset)->gp[i]).estrain[j]=((shells+offset)->gp[i]).estrain[j];
	  ((mshell+offset)->gp[i]).pstrain[j]=((shells+offset)->gp[i]).pstrain[j];
	  ((mshell+offset)->gp[i]).stress[j]=((shells+offset)->gp[i]).stress[j];
	  ((mshell+offset)->gp[i]).backstress[j]=((shells+offset)->gp[i]).backstress[j];
	}
	((mshell+offset)->gp[i]).qn=((shells+offset)->gp[i]).qn;
	((mshell+offset)->gp[i]).qm=((shells+offset)->gp[i]).qm;
	((mshell+offset)->gp[i]).qnm=((shells+offset)->gp[i]).qnm;
	((mshell+offset)->gp[i]).qilyushin=((shells+offset)->gp[i]).qilyushin;
	((mshell+offset)->gp[i]).yinit=((shells+offset)->gp[i]).yinit;
	((mshell+offset)->gp[i]).y=((shells+offset)->gp[i]).y;

	((mshell+offset)->gp[i]).alpha=((shells+offset)->gp[i]).alpha;
	((mshell+offset)->gp[i]).Ee=((shells+offset)->gp[i]).Ee;
	((mshell+offset)->gp[i]).Ep=((shells+offset)->gp[i]).Ep;
  }
  return;
}/*inputmemoryshell*/

/*ELASTO-PLASTICITY*/
double* ilyushin(struct oshell* shell, int ii, double* lambda, double** W)
{
  int i,j;
  int nstress = shell->nstress;
  struct gausspoint* gp;
  double* qstress,* ostress;
  double* f;
  double fn,fm;
  double fn2,fm2,fnfm;
  double E,t,A,I,v;
  double yinit,y;
  double alpha;
  double qn,qm,qnm;

  double c[3],p[3],e[3],det[3];
  double O[2][2][3], Oinv[2][2][3];

  char str[800];

  gp = &(shell->gp[ii]);

  qstress=(double*)malloc(nstress * sizeof(double));
  ostress=(double*)malloc(nstress * sizeof(double));
  f=(double*)malloc(2 * sizeof(double));

  /*          | 1 -1  0 |  */
  /*[Q] = 1/ã2| 1  1  0 |  */
  /*          | 0  0 ã2 |  */

  /*[Q]diag[a,b,c][Q]^T  */
  /*    | a+b a-b   0 |  */
  /*=1/2| a+b a-b   0 |  */
  /*    |   0   0  2c |  */

  /*{ƒÐtry-ƒ¿try}=([I]+2ƒ¢ƒÉ[C][A]){ƒÐ-ƒ¿}=[W]^-1{ƒÐ-ƒ¿} */
  /*{ƒÐ-ƒ¿}=[W]{ƒÐtry-ƒ¿try}                         */
  /*[W]^-1=(diag[Q,Q])[O](diag[Q^T,Q^T])         */
  /*[W]=(diag[Q,Q])[O]^-1(diag[Q^T,Q^T])         */

  E=shell->sect->E;
  t=shell->sect->area;
  v=shell->sect->poi;
  A=E*t;
  I=E*pow(t,3)/12.0;
  yinit=shell->sect->yieldinit;

  fn=shell->sect->fmax[0] /*yinit*t*/;
  fm=shell->sect->fmax[3] /*yinit*pow(t,2)/4.0*/;
  /**/

  fn2=pow(fn,2);
  fm2=pow(fm,2);
  fnfm=fn*fm;


  c[0]=1/(1-v);
  c[1]=1/(1+v);
  c[2]=1/(2.0*(1+v));

  p[0]=0.5;
  p[1]=1.5;
  p[2]=3.0;

  e[0] = (*(lambda+0) + *(lambda+1)) / fn2;
  e[1] = (*(lambda+0) + *(lambda+1)) / fm2;
  e[2] = (*(lambda+0) - *(lambda+1)) /(2.0*sqrt(3.0)*fnfm);


  /*[O]*/
  for (i = 0; i < 3; i++)
  {
	O[0][0][i] = 2.0*A*c[i]*p[i] /*+ 4.0/3.0*Hkin*[p]*/; O[0][1][i] = O[0][0][i];
	O[1][0][i] = 2.0*I*c[i]*p[i] /*+ 4.0/3.0*Hkin*[p]*/; O[1][1][i] = O[1][0][i];
  }

  for (i = 0; i < 3; i++)
  {
	O[0][0][i] = e[0]*O[0][0][i] + 1.0;
	O[0][1][i] = e[2]*O[0][1][i];
	O[1][0][i] = e[2]*O[1][0][i];
	O[1][1][i] = e[1]*O[1][1][i] + 1.0;
  }

  /*Oinv:[O]^-1*/
  for (i = 0; i < 3; i++)
  {
	det[i] = O[0][0][i]*O[1][1][i]-O[0][1][i]*O[1][0][i];

	Oinv[0][0][i] =  O[1][1][i]/det[i];
	Oinv[0][1][i] = -O[0][1][i]/det[i];
	Oinv[1][0][i] = -O[1][0][i]/det[i];
	Oinv[1][1][i] =  O[0][0][i]/det[i];
  }

  //if(shell->loff==20 && ii==0)dbgvct(Oinv[0][0],3,3,"Oinv00");
  //if(shell->loff==20 && ii==0)dbgvct(Oinv[0][1],3,3,"Oinv01");
  //if(shell->loff==20 && ii==0)dbgvct(Oinv[1][0],3,3,"Oinv10");
  //if(shell->loff==20 && ii==0)dbgvct(Oinv[1][1],3,3,"Oinv11");


  /*(diag[Q^T,Q^T]){ƒÐtry-ƒ¿try}*/
  *(qstress+0)=( gp->stress[0] + gp->stress[1] )/sqrt(2.0);
  *(qstress+1)=(-gp->stress[0] + gp->stress[1] )/sqrt(2.0);
  *(qstress+2)=  gp->stress[2];
  *(qstress+3)=( gp->stress[3] + gp->stress[4] )/sqrt(2.0);
  *(qstress+4)=(-gp->stress[3] + gp->stress[4] )/sqrt(2.0);
  *(qstress+5)=  gp->stress[5];

  /*[O]^-1*(diag[Q^T,Q^T]){ƒÐtry-ƒ¿try}*/
  *(ostress+0) = Oinv[0][0][0]**(qstress+0) + Oinv[0][1][0]**(qstress+3);
  *(ostress+1) = Oinv[0][0][1]**(qstress+1) + Oinv[0][1][1]**(qstress+4);
  *(ostress+2) = Oinv[0][0][2]**(qstress+2) + Oinv[0][1][2]**(qstress+5);
  *(ostress+3) = Oinv[1][0][0]**(qstress+0) + Oinv[1][1][0]**(qstress+3);
  *(ostress+4) = Oinv[1][0][1]**(qstress+1) + Oinv[1][1][1]**(qstress+4);
  *(ostress+5) = Oinv[1][0][2]**(qstress+2) + Oinv[1][1][2]**(qstress+5);



  /*f={ƒÐtry-ƒ¿try}^T[W]^T[A][W]{ƒÐtry-ƒ¿try}*/
  qn = (p[0]**(ostress+0)**(ostress+0) + p[1]**(ostress+1)**(ostress+1) + p[2]**(ostress+2)**(ostress+2))/fn2 ;
  qm = (p[0]**(ostress+3)**(ostress+3) + p[1]**(ostress+4)**(ostress+4) + p[2]**(ostress+5)**(ostress+5))/fm2 ;
  qnm= (p[0]**(ostress+0)**(ostress+3) + p[1]**(ostress+1)**(ostress+4) + p[2]**(ostress+2)**(ostress+5))/fnfm;

  *(f+0) = qn+qm+(qnm/sqrt(3.0));
  *(f+1) = qn+qm-(qnm/sqrt(3.0));

  gp->qn = qn;
  gp->qm = qm;
  gp->qnm = qnm;
  gp->qilyushin = *(f+0);

#if 0

  *(f+0) = p[0] * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 + *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + p[1] * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 + *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + p[2] * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 + *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm));
  *(f+1) = p[0] * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 - *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + p[1] * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 - *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + p[2] * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 - *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm));


  qn  = ( pow(gp->stress[0],2)+pow(gp->stress[1],2)-(gp->stress[0])*(gp->stress[1])+3.0*pow(gp->stress[2],2) )/fn2;
  qm  = ( pow(gp->stress[3],2)+pow(gp->stress[4],2)-(gp->stress[3])*(gp->stress[4])+3.0*pow(gp->stress[5],2) )/fm2;
  qnm = (       (gp->stress[0])*(gp->stress[3])
		  +     (gp->stress[1])*(gp->stress[4])
		  -0.5* (gp->stress[0])*(gp->stress[4])
		  -0.5* (gp->stress[1])*(gp->stress[3])
		  +3.0* (gp->stress[2])*(gp->stress[5])  )/fnfm;

#endif

  alpha = gp->alpha + 2.0*( *(lambda+0)*sqrt(*(f+0)) + *(lambda+1)*sqrt(*(f+1)) )/(t*yinit);
  y = yieldstress(shell->sect,alpha,NULL,NULL);
  *(f+0) -= pow(y/yinit,2.0);
  *(f+1) -= pow(y/yinit,2.0);

  gp->yinit = yinit;
  gp->y = y;

  if(W!=NULL)
  {
	  /*(diag[Q,Q])[O]^-1*(diag[Q^T,Q^T]){ƒÐtry-ƒ¿try}*/
	  gp->stress[0]=(*(ostress+0)-*(ostress+1))/sqrt(2.0);
	  gp->stress[1]=(*(ostress+0)+*(ostress+1))/sqrt(2.0);
	  gp->stress[2]= *(ostress+2);
	  gp->stress[3]=(*(ostress+3)-*(ostress+4))/sqrt(2.0);
	  gp->stress[4]=(*(ostress+3)+*(ostress+4))/sqrt(2.0);
	  gp->stress[5]= *(ostress+5);

	  //gp->backstress[0]=(*(ostress+0)-*(ostress+1))/sqrt(2.0);
	  //gp->backstress[1]=(*(ostress+0)+*(ostress+1))/sqrt(2.0);
	  //gp->backstress[2]= *(ostress+2);
	  //gp->backstress[3]=(*(ostress+3)-*(ostress+4))/sqrt(2.0);
	  //gp->backstress[4]=(*(ostress+3)+*(ostress+4))/sqrt(2.0);
	  //gp->backstress[5]= *(ostress+5);

	  gp->alpha = alpha;

	  for (i = 0; i < 2; i++)
	  {
		  for (j = 0; j < 2; j++)
		  {
			  *(*(W+3*i+0)+3*j+0) = 0.5*(Oinv[i][j][0]+Oinv[i][j][1]);
			  *(*(W+3*i+0)+3*j+1) = 0.5*(Oinv[i][j][0]-Oinv[i][j][1]);

			  *(*(W+3*i+1)+3*j+0) = *(*(W+3*i+0)+3*j+1);
			  *(*(W+3*i+1)+3*j+1) = *(*(W+3*i+0)+3*j+0);

			  *(*(W+3*i+2)+3*j+2) = Oinv[i][j][2];
		  }
	  }
	   *(*(W+6)+6) = 1.0;

  }
  free(qstress);
  free(ostress);

  return f;
}

double yieldstress(struct osect* sect, double alpha, double* dy, double* ddy)
{
	char str[256];
	double y;
	double yinit,c,n;
	yinit = sect->yieldinit;
	c = sect->yieldcoefficient;
	n = sect->yieldpower;

	if(n==0.0)
	{
	  sprintf(str,"INVALID MATERIAL : YIELD STRESS FUNCTION");
	  MessageBox(NULL,str,"ERROR",MB_OK);
	}
	else
	{
	  y = yinit + c * pow(alpha,n);
	  if(n==1.0)
	  {
		if(dy!=NULL)*dy = c;
		if(ddy!=NULL)*ddy = 0.0;
	  }
	  else
	  {
		if(dy!=NULL)*dy = c*n*pow(alpha,n-1);
		if(ddy!=NULL)*ddy = c*n*n*pow(alpha,n-2);
	  }
	}

	return y;
}


void assemshellestress(struct oshell* shell, double*** C)
{
  int ii,i,j;
  int nnod,ngp,nstress;
  nnod = shell->nnod;
  ngp = shell->ngp;
  nstress = shell->nstress;

  struct gausspoint* gp;
  double* gpstrain;
  double* gpstress;
  double* lambda;
  double** W;

  double** H, Halpha;
  double** g, galpha;
  double gtotal[nstress];
  double** consistentC;
  double c[3],a[6];




  for(ii = 0; ii < ngp; ii++)/*FOR EACH INTEGRATION POINT*/
  {
	gp = &(shell->gp[ii]);

	gpstrain = (double*)malloc(nstress * sizeof(double));
	for(i = 0; i < nstress; i++)
	{
	  gp->estrain[i] -= gp->pstrain[i];/*TRIAL STRESS*//*pstrain : CONVERGED PLASTIC STRAIN AT LAST LAP.*/
	  *(gpstrain+i) = gp->estrain[i];
	}
	gpstress = matrixvector(*(C+ii),gpstrain,nstress);/*ELASTIC PREDICTOR*/
	for(i = 0; i < nstress; i++)
	{
	  gp->stress[i] = *(gpstress + i);/*TRIAL STRESS*/
	}
	free(gpstress);
	free(gpstrain);

	W = (double**)malloc(nstress * sizeof(double*));
	for (i = 0; i < nstress; i++)
	{
	  *(W + i) = (double*)malloc(nstress * sizeof(double));
	  for (j = 0; j < nstress; j++)
	  {
		*(*(W + i) + j) = 0.0;
	  }
	}

	/*UPDATE STRESS RESULTANT*/
	lambda = returnmapilyushin(shell, ii, W);

	H = matrixmatrix(W,*(C+ii),nstress);

	/*{df/dƒÐ}={d(ƒÐAƒÐ)/dƒÐ}=2[A]{ƒÐ}={g}*/
	double t;
	double fn,fm;
	double y, yinit;
	double dy, ddy;
	yinit = shell->sect->yieldinit;
	t     = shell->sect->area;

	fn=shell->sect->fmax[0]/*yinit*t*/;
	fm=shell->sect->fmax[3]/*yinit*pow(t,2)/4.0*/;

	c[0] = fn*fn;
	c[1] = fm*fm;
	c[2] = 2.0*sqrt(3.0)*fn*fm;

	a[0] = 2.0 * gp->stress[0] - gp->stress[1];
	a[1] = 2.0 * gp->stress[1] - gp->stress[0];
	a[2] = 6.0 * gp->stress[2];
	a[3] = 2.0 * gp->stress[3] - gp->stress[4];
	a[4] = 2.0 * gp->stress[4] - gp->stress[3];
	a[5] = 6.0 * gp->stress[5];

	g  = (double**)malloc(2 * sizeof(double*));
	for (i = 0; i < 2; i++)
	{
	   *(g+i) = (double*)malloc(nstress * sizeof(double));
	}

	*(*(g+0)+0) = a[0]/c[0] + a[3]/c[2];
	*(*(g+0)+1) = a[1]/c[0] + a[4]/c[2];
	*(*(g+0)+2) = a[2]/c[0] + a[5]/c[2];
	*(*(g+0)+3) = a[0]/c[2] + a[3]/c[1];
	*(*(g+0)+4) = a[1]/c[2] + a[4]/c[1];
	*(*(g+0)+5) = a[2]/c[2] + a[5]/c[1];
	*(*(g+0)+6) = 0.0;

	*(*(g+1)+0) =  a[0]/c[0] - a[3]/c[2];
	*(*(g+1)+1) =  a[1]/c[0] - a[4]/c[2];
	*(*(g+1)+2) =  a[2]/c[0] - a[5]/c[2];
	*(*(g+1)+3) = -a[0]/c[2] + a[3]/c[1];
	*(*(g+1)+4) = -a[1]/c[2] + a[4]/c[1];
	*(*(g+1)+5) = -a[2]/c[2] + a[5]/c[1];
	*(*(g+1)+6) =  0.0;

	/*UPDATE HARDENING PARAMS*/
	for (i = 0; i < nstress; i++)
	{
	  gtotal[i] = *(lambda+0)**(*(g+0)+i) + *(lambda+1)**(*(g+1)+i);
	  gp->estrain[i] -= gtotal[i];
	  gp->pstrain[i] += gtotal[i];
	}

	/*df/dq=d/dq(-ƒÐy^2/ƒÐyinit^2)=-2*ƒÐy*ƒÐy'(q)/(ƒÐyinit^2)=beta*/

	y = yieldstress(shell->sect, gp->alpha, &dy, &ddy);
	galpha = -2.0*y*dy/pow(yinit,2.0);
	Halpha = 1.0/(t*dy-(*(lambda+0)+*(lambda+1))*2.0*(dy*dy+y*ddy)/pow(yinit,2.0));

	/*UPDATE CONSISTENT STIFFNESS*/
	consistentC = consistentCilyushin(shell, ii, H, Halpha, g, galpha, nstress);
	for (i = 0; i < nstress; i++)
	{
	  for (j = 0; j < nstress; j++)
	  {
		*(*(*(C+ii) + i) + j) = *(*(consistentC + i) + j);
	  }
	}

	freematrix(consistentC,nstress);
	freematrix(W,nstress);
	free(lambda);
  }
  return;
}

#if 0
double* returnmapilyushin(struct oshell* shell, int ii, double** W)
{
   int i,j;
   //struct gausspoint* gp;
   double* lambda,* lambdaeps,* dlambda, *lastlambda;
   double* f,* feps;
   double** dfdl;
   double det;
   double residual = 1.0;
   double tolerance = 1.0e-8;
   double eps = 1.0;
   int FLAG;
   char str[800];



   //gp = &(shell->gp[ii]);
   /*FLAG = 0 ; BOTH SURFACES ARE ACTIVE*/
   /*FLAG = 1 ; SURFACE 1 IS ACTIVE     */
   /*FLAG = 2 ; SURFACE 2 IS ACTIVE     */

   lastlambda = (double*)malloc(2 * sizeof(double));
   lambda = (double*)malloc(2 * sizeof(double));
   dlambda = (double*)malloc(2 * sizeof(double));
   lambdaeps = (double*)malloc(2 * sizeof(double));

   dfdl= (double**)malloc(2 * sizeof(double*));
   for(i=0;i<2;i++)
   {
	 *(dfdl+i)= (double*)malloc(2 * sizeof(double));
   }

   for(i=0;i<2;i++)
   {
	 *(lambda+i) = 0.0;
   }
   f=ilyushin(shell, ii, lambda, W);/*ILYUSHIN'S YIELD FUNCTION.{ƒÐ-ƒ¿}^T[A]{ƒÐ-ƒ¿}*/


   if( *(f+0)>0 || *(f+1)>0 )/*IF BOTH YIELD SURFACES ARE ACTIVE.*/
   {
	   //sprintf(str,"YIELD DETECTED SHELL %d INTEGRATION POINT %d\n",shell->loff,ii);
	   //dbgstr(str);

	   //sprintf(str,"%f %f %f %f %f",(shell->gp[ii]).qn,(shell->gp[ii]).qm,(shell->gp[ii]).qnm,(shell->gp[ii]).yinit,(shell->gp[ii]).y);
	   //dbgstr(str);




	   /*RETURN-MAPPING PROCEDURE USING NEWTON-RAPTHON METHOD FOR ILYUSHIN'S YIELD CONDITION.*/
	   FLAG = 0;
	   while(residual>tolerance)/*BOTH YIELD SURFACES ARE ACTIVE.*/
	   {
		 /*MEMORY*/
		 *(lastlambda + 0) = *(lambda + 0);
		 *(lastlambda + 1) = *(lambda + 1);


		 /*SENSITIVITY*/
		 for(i=0;i<2;i++)
		 {
		   for(j=0;j<2;j++)
		   {
			 if(i==j) *(lambdaeps+j)=*(lambda+j)+eps;
			 else     *(lambdaeps+j)=*(lambda+j);
		   }
		   feps=ilyushin(shell, ii, lambdaeps, NULL);
		   for(j=0;j<2;j++)
		   {
			 *(*(dfdl+j)+i)=(*(feps+j)-*(f+j))/eps;
		   }
		   free(feps);
		 }

		 /*UPDATE*/
		 det=*(*(dfdl+0)+0)**(*(dfdl+1)+1)-*(*(dfdl+0)+1)**(*(dfdl+1)+0);

		 *(dlambda + 0) = -( *(*(dfdl+1)+1)**(f+0)-*(*(dfdl+0)+1)**(f+1))/det;
		 *(dlambda + 1) = -(-*(*(dfdl+1)+0)**(f+0)+*(*(dfdl+0)+0)**(f+1))/det;

		 *(lambda + 0) += *(dlambda + 0);
		 *(lambda + 1) += *(dlambda + 1);

		 /*NEW JUDGE*/
		 free(f);
		 f=ilyushin(shell, ii, lambda, NULL);

		 //sprintf(str,"%f %f %d",*(f+0),*(f+1),FLAG);
		 //dbgstr(str);



		 if(*(f+0) < 0)
		 {
		   /*SECOND YIELD SURFACE IS ACTIVE*/
		   FLAG = 2;
		   *(lambda + 0) = 0.0;
		   *(lambda + 1) = *(lastlambda + 1);
		   residual = *(f+1);
		   break;
		 }
		 if(*(f+1) < 0)
		 {
		   /*FIRST YIELD SURFACE IS ACTIVE*/
		   FLAG = 1;
		   *(lambda + 0) = *(lastlambda + 0);
		   *(lambda + 1) = 0.0;
		   residual = *(f+0);
		   break;
		 }
		 else
		 {
		   residual=vectorlength(f,2);
		 }
	   }

	   if(FLAG == 1)
	   {
		 while(residual>tolerance)/*FIRST YIELD SURFACE IS ACTIVE.*//*f1=0 && f2<0 && lambda1>0 && lambda2=0*/
		 {
		   /*MEMORY*/
		   *(lastlambda + 0) = *(lambda + 0);
		   *(lastlambda + 1) = *(lambda + 1);

		   /*SENSITIVITY*/
		   *(lambdaeps+0)=*(lambda+0)+eps;
		   *(lambdaeps+1)=0.0;
		   feps=ilyushin(shell, ii, lambdaeps, NULL);

		   /*UPDATE*/
		   *(dlambda+0) = -*(f+0)*eps / (*(feps+0)-*(f+0));
		   *(dlambda+1) = 0.0;
		   free(feps);

		   *(lambda+0) += *(dlambda+0);
		   *(lambda+1) += *(dlambda+1);

		   /*NEW JUDGE*/
		   free(f);
		   f=ilyushin(shell, ii, lambda, NULL);

		   sprintf(str,"%f %f %d",*(f+0),*(f+1),FLAG);
		   dbgstr(str);



		   residual = *(f+0);
		 }
	   }
	   if(FLAG == 2)
	   {
		 while(residual>tolerance)/*YIELD SURFACE 2 IS ACTIVE.*//*f1<0 && f2=0 && lambda1=0 && lambda2>0*/
		 {
		   /*MEMORY*/
		   *(lastlambda + 0) = *(lambda + 0);
		   *(lastlambda + 1) = *(lambda + 1);

		   /*SENSITIVITY*/
		   *(lambdaeps+0)=0.0;
		   *(lambdaeps+1)=*(lambda+1)+eps;
		   feps=ilyushin(shell, ii, lambdaeps, NULL);

		   /*UPDATE*/
		   *(dlambda+0) = 0.0;
		   *(dlambda+1) = -*(f+1)*eps / (*(feps+1)-*(f+1));
		   free(feps);

		   *(lambda+0) += *(dlambda+0);
		   *(lambda+1) += *(dlambda+1);

		   /*NEW JUDGE*/
		   free(f);
		   f=ilyushin(shell, ii, lambda, NULL);


		   sprintf(str,"%f %f %d",*(f+0),*(f+1),FLAG);
		   dbgstr(str);


		   residual = *(f+1);
		 }
	   }


	   free(lastlambda);
	   free(lambdaeps);
	   free(dlambda);
	   freematrix(dfdl,2);
   }

   sprintf(str,"CONVERGED %f %f %f %f %f\n\n",(shell->gp[ii]).qn,(shell->gp[ii]).qm,(shell->gp[ii]).qnm,(shell->gp[ii]).yinit,(shell->gp[ii]).y);
   dbgstr(str);



  free(f);

   return lambda;
}
#endif





double* returnmapilyushin(struct oshell* shell, int ii, double** W)
{
   int i,j;
   //struct gausspoint* gp;
   double* lambda,* lambda01,* lambda0,* lambda1,* lambdaeps,* dlambda;
   double* f,* feps,* finit,* f01,* f0,* f1;
   double** dfdl;
   double det;
   double residual;
   double tolerance = 1.0e-8;
   double eps = 1.0e-2;
   char str[800];
   int iteration;



   //gp = &(shell->gp[ii]);
   /*FLAG = 0 ; BOTH SURFACES ARE ACTIVE*/
   /*FLAG = 1 ; SURFACE 1 IS ACTIVE     */
   /*FLAG = 2 ; SURFACE 2 IS ACTIVE     */

   lambda = (double*)malloc(2 * sizeof(double));
   lambda01 = (double*)malloc(2 * sizeof(double));
   lambda0  = (double*)malloc(2 * sizeof(double));
   lambda1  = (double*)malloc(2 * sizeof(double));

   finit = (double*)malloc(2 * sizeof(double));
   f01 = (double*)malloc(2 * sizeof(double));
   f0  = (double*)malloc(2 * sizeof(double));
   f1  = (double*)malloc(2 * sizeof(double));

   dlambda = (double*)malloc(2 * sizeof(double));
   lambdaeps = (double*)malloc(2 * sizeof(double));

   dfdl= (double**)malloc(2 * sizeof(double*));
   for(i=0;i<2;i++)
   {
	 *(dfdl+i)= (double*)malloc(2 * sizeof(double));
   }

   *(lambda + 0) = 0.0;
   *(lambda + 1) = 0.0;
   f=ilyushin(shell, ii, lambda, W);/*ILYUSHIN'S YIELD FUNCTION.{ƒÐ-ƒ¿}^T[A]{ƒÐ-ƒ¿}*/
   *(finit + 0) = *(f + 0);
   *(finit + 1) = *(f + 1);


   if(*(f+0)>0 || *(f+1)>0)/*RETURN-MAPPING PROCEDURE USING NEWTON-RAPTHON METHOD FOR ILYUSHIN'S YIELD CONDITION.*/
   {
	   /*
	   sprintf(str,"YIELD DETECTED SHELL %d INTEGRATION POINT %d\n",shell->loff,ii);
	   dbgstr(str);

	   sprintf(str,"%f %f %f %f %f",(shell->gp[ii]).qn,(shell->gp[ii]).qm,(shell->gp[ii]).qnm,(shell->gp[ii]).yinit,(shell->gp[ii]).y);
	   dbgstr(str);
       */


	   /*BOTH YIELD SURFACES ARE ACTIVE.*/
	   *(lambda + 0) = 0.0;
	   *(lambda + 1) = 0.0;

	   *(f + 0) = *(finit + 0);
	   *(f + 1) = *(finit + 1);

	   residual=1.0;
	   iteration=1;
	   while(residual>tolerance && iteration<10)
	   {
		 /*SENSITIVITY*/
		 for(i=0;i<2;i++)
		 {
		   for(j=0;j<2;j++)
		   {
			 if(i==j) *(lambdaeps+j)=*(lambda+j)+eps;
			 else     *(lambdaeps+j)=*(lambda+j);
		   }
		   feps=ilyushin(shell, ii, lambdaeps, NULL);
		   for(j=0;j<2;j++)
		   {
			 *(*(dfdl+j)+i)=(*(feps+j)-*(f+j))/eps;
		   }
		   free(feps);
		 }

		 /*UPDATE*/
		 det=*(*(dfdl+0)+0)**(*(dfdl+1)+1)-*(*(dfdl+0)+1)**(*(dfdl+1)+0);

		 *(dlambda + 0) = -( *(*(dfdl+1)+1)**(f+0)-*(*(dfdl+0)+1)**(f+1))/det;
		 *(dlambda + 1) = -(-*(*(dfdl+1)+0)**(f+0)+*(*(dfdl+0)+0)**(f+1))/det;

		 *(lambda + 0) += *(dlambda + 0);
		 *(lambda + 1) += *(dlambda + 1);

		 /*NEW JUDGE*/
		 free(f);
		 f=ilyushin(shell, ii, lambda, NULL);

		 //sprintf(str,"%f %f %f %f %d",*(f+0),*(f+1),*(lambda+0),*(lambda+1),0);
		 //dbgstr(str);

		 residual=vectorlength(f,2);
		 iteration++;
	   }
	   *(lambda01 + 0) = *(lambda + 0);
	   *(lambda01 + 1) = *(lambda + 1);

	   *(f01 + 0) = *(f + 0);
	   *(f01 + 1) = *(f + 1);

	   /*FIRST YIELD SURFACE IS ACTIVE.*//*f1=0 && f2<0 && lambda1>0 && lambda2=0*/
	   *(lambda + 0) = 0.0;
	   *(lambda + 1) = 0.0;

	   *(f + 0) = *(finit + 0);
	   *(f + 1) = *(finit + 1);

	   residual=1.0;
	   iteration=1;
	   while(residual>tolerance && iteration<10)
	   {
		 /*SENSITIVITY*/
		 *(lambdaeps+0)=*(lambda+0)+eps;
		 feps=ilyushin(shell, ii, lambdaeps, NULL);

		 /*UPDATE*/
		 *(dlambda+0) = -*(f+0)*eps / (*(feps+0)-*(f+0));
		 free(feps);

		 *(lambda+0) += *(dlambda+0);

		 /*NEW JUDGE*/
		 free(f);
		 f=ilyushin(shell, ii, lambda, NULL);

		 //sprintf(str,"%f %f %f %f %d",*(f+0),*(f+1),*(lambda+0),*(lambda+1),1);
		 //dbgstr(str);

		 residual = abs(*(f+0));
		 iteration++;
	   }
	   *(lambda0 + 0) = *(lambda + 0);
	   *(lambda0 + 1) = *(lambda + 1);

	   *(f0 + 0) = *(f + 0);
	   *(f0 + 1) = *(f + 1);

	   /*YIELD SURFACE 2 IS ACTIVE.*//*f1<0 && f2=0 && lambda1=0 && lambda2>0*/
	   *(lambda + 0) = 0.0;
	   *(lambda + 1) = 0.0;

	   *(f + 0) = *(finit + 0);
	   *(f + 1) = *(finit + 1);

	   residual=1.0;
	   iteration=1;
	   while(residual>tolerance && iteration<10)
	   {
		 /*SENSITIVITY*/
		 *(lambdaeps+1)=*(lambda+1)+eps;
		 feps=ilyushin(shell, ii, lambdaeps, NULL);

		 /*UPDATE*/
		 *(dlambda+1) = -*(f+1)*eps / (*(feps+1)-*(f+1));
		 free(feps);

		 *(lambda+1) += *(dlambda+1);

		 /*NEW JUDGE*/
		 free(f);
		 f=ilyushin(shell, ii, lambda, NULL);

		 //sprintf(str,"%f %f %f %f %d",*(f+0),*(f+1),*(lambda+0),*(lambda+1),2);
		 //dbgstr(str);

		 residual = abs(*(f+1));
		 iteration++;
	   }
	   *(lambda1 + 0) = *(lambda + 0);
	   *(lambda1 + 1) = *(lambda + 1);

	   *(f1 + 0) = *(f + 0);
	   *(f1 + 1) = *(f + 1);



	   if( *(lambda01 + 0) > 0.0 && *(lambda01 + 1) > 0.0 && *(f01 + 0) <  tolerance && *(f01 + 1) <  tolerance)
	   {
		 *(lambda + 0) = *(lambda01 + 0);
		 *(lambda + 1) = *(lambda01 + 1);
	   }
	   if( *(lambda0 + 0) > 0.0 && abs(*(f0 + 0)) < tolerance && *(f0 + 1) < 0.0)
	   {
		 *(lambda + 0) = *(lambda0 + 0);
		 *(lambda + 1) = *(lambda0 + 1);
	   }
	   if( *(lambda1 + 1) > 0.0 && abs(*(f1 + 1)) <  tolerance && *(f1 + 0) < 0.0)
	   {
		 *(lambda + 0) = *(lambda1 + 0);
		 *(lambda + 1) = *(lambda1 + 1);
	   }

	  free(f);
	  f=ilyushin(shell, ii, lambda, W);

	  //sprintf(str,"CONVERGED %f %f %f %f %f %f %f\n\n",*(lambda + 0),*(lambda + 1),(shell->gp[ii]).qn,(shell->gp[ii]).qm,(shell->gp[ii]).qnm,(shell->gp[ii]).yinit,(shell->gp[ii]).y);
	  //dbgstr(str);
   }

	  (shell->gp[ii]).f[0]=*(f+0);
	  (shell->gp[ii]).f[1]=*(f+1);



   free(lambdaeps);
   free(dlambda);
   freematrix(dfdl,2);

   free(f01);
   free(f0);
   free(f1);
   free(f);
   free(finit);
   free(lambda01);
   free(lambda0);
   free(lambda1);


   return lambda;
}



double** consistentCilyushin(struct oshell* shell, int ii, double** H, double Halpha, double** g, double galpha, int nstress)
{
  /*FOR ONE GAUSS INTEGRATION POINT*/
  int i,j;
  double** consistentC;
  double* g0,* g1,* N0,* N1;
  double M00,M01,M10,M11;
  double det;
  double beta;
  double tolerance = 1.0e-8;

  consistentC = (double**)malloc(nstress * sizeof(double*));
  for (i = 0; i < nstress; i++)
  {
	*(consistentC + i) = (double*)malloc(nstress * sizeof(double));
	for (j = 0; j < nstress; j++)
	{
		*(*(consistentC+i)+j) = *(*(H+i)+j);
	}
  }

  beta=Halpha*pow(galpha,2.0);

  if(abs((shell->gp[ii]).f[0])<tolerance && abs((shell->gp[ii]).f[1])<tolerance)
  {
	N0=matrixvector(H,*(g+0),nstress);
	N1=matrixvector(H,*(g+1),nstress);

	M00=dotproduct(*(g+0),N0,nstress);
	M01=dotproduct(*(g+0),N1,nstress);
	M10=M01;
	M11=dotproduct(*(g+1),N1,nstress);

	M00+=beta;
	M01+=beta;
	M10+=beta;
	M11+=beta;
	det = M00*M11-M01*M10;

	for (i = 0; i < nstress; i++)
	{
	  for (j = 0; j < nstress; j++)
	  {
		*(*(consistentC + i) + j) -= ( M11**(N0+i)**(N0+j) - M01**(N0+i)**(N1+j) - M10**(N1+i)**(N0+j) + M00**(N1+i)**(N1+j) )/det;
	  }
	}
	free(N0);
	free(N1);
  }
  else if(abs((shell->gp[ii]).f[0])<tolerance)
  {
	N0=matrixvector(H,*(g+0),nstress);
	M00=dotproduct(*(g+0),N0,nstress);
	M00+=beta;

	for (i = 0; i < nstress; i++)
	{
	  for (j = 0; j < nstress; j++)
	  {
		*(*(consistentC + i) + j) -= *(N0+i)**(N0+j)/M00;
	  }
	}
	free(N0);
  }
  else if(abs((shell->gp[ii]).f[1])<tolerance)
  {
	N1=matrixvector(H,*(g+1),nstress);
	M11=dotproduct(*(g+1),N1,nstress);
	M11+=beta;

	for (i = 0; i < nstress; i++)
	{
	  for (j = 0; j < nstress; j++)
	  {
		*(*(consistentC + i) + j) -= *(N1+i)**(N1+j)/M11;
	  }
	}
	free(N1);
  }

  return consistentC;
}


void assemshellestrain(struct oshell* shell, double*** B, double* edisp)
{
  int ii,i;
  int nnod,ngp,nstress;
  double* gpstrain;

  nnod = shell->nnod;
  ngp = shell->ngp;
  nstress = shell->nstress;
  for(ii = 0; ii < ngp; ii++)/*FOR EACH INTEGRATION POINT*/
  {
	gpstrain = matrixvectorIII(*(B+ii),edisp,nstress,6*nnod); /*INCREMENTAL STRAIN d{ƒÃx ƒÃy ƒÃxy ƒÁx ƒÁy ƒÁxy} FOR NODE ii*/

	for(i = 0; i < nstress; i++)
	{
	  (shell->gp[ii]).estrain[i] = *(gpstrain+i);
	}
	free(gpstrain);
  }
  return;
}

double* assemshelleinternal(struct oshell* shell, double*** B)
{
  int ii,i,j;
  int nnod,ngp,nstress;
  double* einternal,*gpstress,*gpinternal;
  double** Bt;

  nnod = shell->nnod;
  ngp = shell->ngp;
  nstress = shell->nstress;

  einternal = (double*)malloc(6 * nnod * sizeof(double));
  for(i = 0; i < nnod; i++)
  {
	for(j = 0; j < 6; j++)
	{
	   *(einternal+6*i+j) = 0.0/*shell->stress[i][j]*/;
	}
  }
  gpstress = (double*)malloc(nstress * sizeof(double));

  for(ii = 0; ii < ngp; ii++)
  {
	for(i = 0; i < nstress; i++)
	{
	  *(gpstress+i) = (shell->gp[ii]).stress[i];
	}
	Bt = matrixtransposeIII(*(B+ii),nstress,6*nnod);
	gpinternal = matrixvectorIII(Bt,gpstress,6*nnod,nstress);
	for(i = 0; i < 6*nnod; i++)
	{
	  *(einternal+i) += *(gpinternal+i)*(shell->w[ii])*(shell->area);
	}
	freematrix(Bt,6*nnod);
	free(gpinternal);
  }
  free(gpstress);
  for(i = 0; i < nnod; i++)
  {
	for(j = 0; j < 6; j++)
	{
	  shell->stress[i][j] = *(einternal+6*i+j);
	}
  }
  return einternal;
}







