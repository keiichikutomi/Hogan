
/*FUNCTIONS FOR SHELL*/
double** planeC(struct oshell shell);
double** shellC(struct oshell shell);


double*** shellB(struct oshell shell);
double** shelldrccos(struct oshell shell);
double shellarea(struct oshell shell);
double* extractshelldisplacement(struct oshell shell, double* ddisp);
double** shelllocalcoord(struct oshell shell);


double** assemshellemtx(struct oshell shell);/*ELASTIC STIFFNESS MATRIX*/

double** assemshellpmtx(struct oshell shell, double*** C, double*** B);/*ELASTO-PLASTIC STIFFNESS MATRIX*/

double** assemshellmmtx(struct oshell shell);/*MASS MATRIX*/

double* assemshellpvct(struct oshell shell);/*ELEMENT LOAD*/

double shellvolume(struct oshell shell, double** drccos);/*VOLUME*/
void assemshellvolume(struct oshell* shells, int nshell, double* ddisp, double* volume);


/*INPUT & OUTPUT SHELL DATA*/
void inputshell(struct oshell *shells,struct memoryshell *mshell,int offset,struct oshell *shell);
void outputshell(struct oshell *shells,int offset,struct oshell *shell);
void initialshell(struct oshell *shells,struct memoryshell *mshell,int nshell);/*copy data from shells to mshell*/
/*FOR PLASTIC (ILYUSHIN'S STRESS RESULTANT)*/

double yieldstress(struct osect* sect, double alpha, double* dy, double* ddy);

double*** shellCconsistent(struct oshell shell);/*UPDATE STIFFNESS MATRIX*/

double* ilyushin(struct oshell* shell, int ii, double* lambda, int UPDATEFLAG);/*YIELD FUNCTION*/
void returnmapilyushin(struct oshell* shell, int ii);/*RETURN-MAPPING*/
double** shellCconsistentilyushin(struct oshell shell, int ii);/*UPDATE STIFFNESS MATRIX*/

double mises(struct oshell* shell, int ii, int jj, double* lambda, int UPDATEFLAG);/*YIELD FUNCTION*/
void returnmapmises(struct oshell* shell, int ii, int jj);/*RETURN-MAPPING*/
double** shellCconsistentmises(struct oshell shell, int ii);/*UPDATE STIFFNESS MATRIX*/

void assemshellestrain(struct oshell* shell, double*** B, double* edisp);
void assemshellestress(struct oshell* shell/*, double*** C*/);
double* assemshelleinternal(struct oshell* shell, double*** B);



extern void dbgstr(const char* str);
extern void dbgvct(double* vct, int size, int linesize, const char* str);
extern void dbgmtx(double** mtx, int rows, int cols, const char* str);


double** planeC(struct oshell shell)
{
	int i,j,ii;
	double** C;
	double E = shell.sect->E;
	double t = shell.sect->area;
	double poi = shell.sect->poi;

	//double prate = shell.prate;
	//double brate = shell.brate;

	C = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(C+i) = (double*)malloc(3 * sizeof(double));
		for (j = 0; j < 3; j++)
		{
			*(*(C+i)+j) = 0.0;
		}
	}

	*(*(C+0)+0)=E/(1-poi*poi);
	*(*(C+0)+1)=poi**(*(C+0)+0);
	*(*(C+1)+0)=*(*(C+0)+1);
	*(*(C+1)+1)=*(*(C+0)+0);
	*(*(C+2)+2)=0.5*E/(1+poi);

	return C;
}

double** shellC(struct oshell shell)
{
	int i,j,ii;
	int nstress = shell.nstress;
	double** Cplane;
	double** C;
	double E = shell.sect->E;
	double t = shell.sect->area;
	double poi = shell.sect->poi;

	double prate = shell.prate;
	double brate = shell.brate;


	Cplane = planeC(shell);
	C = (double**)malloc(nstress * sizeof(double*));
	for (i = 0; i < nstress; i++)
	{
		*(C + i) = (double*)malloc(nstress * sizeof(double));
		for (j = 0; j < nstress; j++)
		{
			*(*(C + i) + j) = 0.0;
		}
	}

	/*FOR RESULTANT STRESS*/
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			*(*(C+i)+j)     = *(*(Cplane+i)+j)*t*prate;
			*(*(C+3+i)+3+j) = *(*(Cplane+i)+j)*pow(t,3)/12.0*brate;
		}
	}

	/*
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
	*/
	freematrix(Cplane,3);
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

double** shelllocalcoord(struct oshell shell)
{
  int i,j;
  double** exy,** drccos;

  exy=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)/*NODE 0-0, 0-1, 0-2*/
  {
	*(exy+i)=(double *)malloc(2*sizeof(double));
  }

  drccos = shelldrccos(shell);

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
  return exy;
}

double*** shellB(struct oshell shell)
{
  int i,j,k,ii;
  int nnod = shell.nnod;
  int ngp = shell.ngp;
  int nstress = shell.nstress;
  double **exy;
  double area;
  double** L;
  double *a,*b,*c;
  double *Lx,*Ly;
  double *len,*aa,*bb,*cc,*dd,*ee;
  double alphab;
  double beta[10];
  double** Q1,**Q2,**Q3,**Q;
  double** Th,** Tn,**TnQ;
  double** Bmb,**Bmh,** Bb;
  double*** B;

  /*TRIANGULAR COORDINATION OF INTEGRATION POINTS*/
  L=(double **)malloc(shell.ngp*sizeof(double *));
  for(i=0;i<shell.ngp;i++)
  {
	*(L+i)=(double *)malloc(3*sizeof(double));
  }

  a=(double *)malloc(3*sizeof(double));
  b=(double *)malloc(3*sizeof(double));
  c=(double *)malloc(3*sizeof(double));
  Lx=(double *)malloc(3*sizeof(double));
  Ly=(double *)malloc(3*sizeof(double));

  len=(double *)malloc(3*sizeof(double));
  aa=(double *)malloc(3*sizeof(double));
  bb=(double *)malloc(3*sizeof(double));
  cc=(double *)malloc(3*sizeof(double));
  dd=(double *)malloc(3*sizeof(double));
  ee=(double *)malloc(3*sizeof(double));

  Q1=(double **)malloc(3*sizeof(double *));
  Q2=(double **)malloc(3*sizeof(double *));
  Q3=(double **)malloc(3*sizeof(double *));

  Q=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Q1+i)=(double *)malloc(3*sizeof(double));
	*(Q2+i)=(double *)malloc(3*sizeof(double));
	*(Q3+i)=(double *)malloc(3*sizeof(double));

	*(Q+i)=(double *)malloc(3*sizeof(double));
  }

  Th=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Th+i)=(double *)malloc(9*sizeof(double));
  }
  Tn=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Tn+i)=(double *)malloc(3*sizeof(double));
  }


  Bmb=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bmb+i)=(double *)malloc(9*sizeof(double));
  }
  Bb=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bb+i)=(double *)malloc(9*sizeof(double));
  }

  B=(double ***)malloc(shell.ngp*sizeof(double **));
  for(ii=0;ii<shell.ngp;ii++)
  {
	*(B+ii)=(double **)malloc(nstress*sizeof(double*));
	for(i=0;i<nstress;i++)
	{
	  *(*(B+ii)+i)=(double *)malloc(6*nnod*sizeof(double));
	  for(j=0;j<6*nnod;j++)
	  {
		*(*(*(B+ii)+i)+j) = 0.0;
	  }
	}
  }

  exy = shelllocalcoord(shell);
  area = shell.area;

  /*AREA COORD Li=(ai+bix+ciy)/(2*area)*/
  if(nnod==3 && ngp==7)
  {
	  /*3 NODES 7 GAUSS POINTS*/
	  *(*(L+0)+0)=1.0/3.0;*(*(L+0)+1)=1.0/3.0;*(*(L+0)+2)=1.0/3.0;

	  *(*(L+1)+0)=1.0/2.0;*(*(L+1)+1)=1.0/2.0;*(*(L+1)+2)=0.0;
	  *(*(L+2)+0)=0.0;*(*(L+2)+1)=1.0/2.0;*(*(L+2)+2)=1.0/2.0;
	  *(*(L+3)+0)=1.0/2.0;*(*(L+3)+1)=0.0;*(*(L+3)+2)=1.0/2.0;

	  *(*(L+4)+0)=1.0;*(*(L+4)+1)=0.0;*(*(L+4)+2)=0.0;
	  *(*(L+5)+0)=0.0;*(*(L+5)+1)=1.0;*(*(L+5)+2)=0.0;
	  *(*(L+6)+0)=0.0;*(*(L+6)+1)=0.0;*(*(L+6)+2)=1.0;
  }
  if(nnod==3 && ngp==3)
  {
	  /*3 NODES 3 GAUSS POINTS*/
	  *(*(L+0)+0)=1.0/6.0;*(*(L+0)+1)=1.0/6.0;*(*(L+0)+2)=2.0/3.0;
	  *(*(L+1)+0)=1.0/6.0;*(*(L+1)+1)=2.0/3.0;*(*(L+1)+2)=1.0/6.0;
	  *(*(L+2)+0)=2.0/3.0;*(*(L+2)+1)=1.0/6.0;*(*(L+2)+2)=1.0/6.0;
  }
  if(nnod==3 && ngp==4)
  {
	  /*3 NODES 4 GAUSS POINTS*/
	  *(*(L+0)+0)=1.0/3.0;*(*(L+0)+1)=1.0/3.0;*(*(L+0)+2)=1.0/3.0;
	  *(*(L+1)+0)=1.0/6.0;*(*(L+1)+1)=1.0/6.0;*(*(L+1)+2)=2.0/3.0;
	  *(*(L+2)+0)=1.0/6.0;*(*(L+2)+1)=2.0/3.0;*(*(L+2)+2)=1.0/6.0;
	  *(*(L+3)+0)=2.0/3.0;*(*(L+3)+1)=1.0/6.0;*(*(L+3)+2)=1.0/6.0;
  }

  /*
  6             |\
  | \           | \
  |  \          |  \
  |   \         |   \
  3    2        |2  3\
  |     \       |     \
  |   0  \      |   0  \
  |       \     |    1  \
  4----1----5   ----------
  */

  for(i=0;i<3;i++)
  {
	  j=(i+1)%3;
	  k=(i+2)%3;

	  *(a+i)=*(*(exy+j)+0)**(*(exy+k)+1)-*(*(exy+k)+0)**(*(exy+j)+1);/*xjyk-xkyj*/
	  *(b+i)=*(*(exy+j)+1)-*(*(exy+k)+1);/*yj-yk*/
	  *(c+i)=*(*(exy+k)+0)-*(*(exy+j)+0);/*xk-xj*/

	  *(Lx+i)=0.5**(b+i)/area;/*dLi/dx=bi/(2*area)*/
	  *(Ly+i)=0.5**(c+i)/area;/*dLi/dy=ci/(2*area)*/

	  *(len+i)= *(b+i)**(b+i) + *(c+i)**(c+i);
	  *(aa+i)= *(c+i)                                / *(len+i);
	  *(bb+i)=-0.75**(b+i)**(c+i)                    / *(len+i);
	  *(cc+i)=(0.25**(c+i)**(c+i)-0.5**(b+i)**(b+i)) / *(len+i);
	  *(dd+i)=-*(b+i)                                / *(len+i);
	  *(ee+i)=(0.25**(b+i)**(b+i)-0.5**(c+i)**(c+i)) / *(len+i);
  }


  /*FREE PARAMETERS OF ANDES OPTIMAL TRIANGLE*/

  double poi = shell.sect->poi;
  alphab = 1.5;
  if(poi!=0.5)
  {
	beta[0] = 0.5*(1-4*poi*poi);
  }
  else
  {
	beta[0] = 0.01;//FOR STABILITY
  }
  beta[1] = 1.0; beta[2] = 2.0; beta[3] = 1.0;
  beta[4] = 0.0; beta[5] = 1.0; beta[6] =-1.0;
  beta[7] =-1.0; beta[8] =-1.0; beta[9] =-2.0;


  *(*(Q1+0)+0) = 2.0*area*beta[1]/(3.0**(len+2));
  *(*(Q1+0)+1) = 2.0*area*beta[2]/(3.0**(len+2));
  *(*(Q1+0)+2) = 2.0*area*beta[3]/(3.0**(len+2));
  *(*(Q1+1)+0) = 2.0*area*beta[4]/(3.0**(len+0));
  *(*(Q1+1)+1) = 2.0*area*beta[5]/(3.0**(len+0));
  *(*(Q1+1)+2) = 2.0*area*beta[6]/(3.0**(len+0));
  *(*(Q1+2)+0) = 2.0*area*beta[7]/(3.0**(len+1));
  *(*(Q1+2)+1) = 2.0*area*beta[8]/(3.0**(len+1));
  *(*(Q1+2)+2) = 2.0*area*beta[9]/(3.0**(len+1));

  *(*(Q2+0)+0) = 2.0*area*beta[9]/(3.0**(len+2));
  *(*(Q2+0)+1) = 2.0*area*beta[7]/(3.0**(len+2));
  *(*(Q2+0)+2) = 2.0*area*beta[8]/(3.0**(len+2));
  *(*(Q2+1)+0) = 2.0*area*beta[3]/(3.0**(len+0));
  *(*(Q2+1)+1) = 2.0*area*beta[1]/(3.0**(len+0));
  *(*(Q2+1)+2) = 2.0*area*beta[2]/(3.0**(len+0));
  *(*(Q2+2)+0) = 2.0*area*beta[6]/(3.0**(len+1));
  *(*(Q2+2)+1) = 2.0*area*beta[4]/(3.0**(len+1));
  *(*(Q2+2)+2) = 2.0*area*beta[5]/(3.0**(len+1));

  *(*(Q3+0)+0) = 2.0*area*beta[5]/(3.0**(len+2));
  *(*(Q3+0)+1) = 2.0*area*beta[6]/(3.0**(len+2));
  *(*(Q3+0)+2) = 2.0*area*beta[4]/(3.0**(len+2));
  *(*(Q3+1)+0) = 2.0*area*beta[8]/(3.0**(len+0));
  *(*(Q3+1)+1) = 2.0*area*beta[9]/(3.0**(len+0));
  *(*(Q3+1)+2) = 2.0*area*beta[7]/(3.0**(len+0));
  *(*(Q3+2)+0) = 2.0*area*beta[2]/(3.0**(len+1));
  *(*(Q3+2)+1) = 2.0*area*beta[3]/(3.0**(len+1));
  *(*(Q3+2)+2) = 2.0*area*beta[1]/(3.0**(len+1));

  /*TRANSFORMATION MATRIX TO EXTRACT HIERARCHICAL ROTATION*/
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Th+i)+3*j+0) = 0.5**(Ly+j);
	  *(*(Th+i)+3*j+1) =-0.5**(Lx+j);
	  if(i==j)
	  {
		*(*(Th+i)+3*j+2) = 1.0;
	  }
	  else
	  {
		*(*(Th+i)+3*j+2) = 0.0;
	  }
	}
  }
  /*STRAIN TRANSFORMATION MATRIX FROM NATURAL COORD TO CARTESIAN COORD*/
  for(i=0;i<3;i++)
  {
	j=(i+1)%3;
	k=(i+2)%3;

	*(*(Tn+0)+i) = -*(Lx+i)**(Lx+j)**(len+k);
	*(*(Tn+1)+i) = -*(Ly+i)**(Ly+j)**(len+k);
	*(*(Tn+2)+i) = -(*(Lx+i)**(Ly+j)+*(Ly+i)**(Lx+j))**(len+k);
  }


  for(i=0;i<3;i++)
  {
	j=(i+1)%3;
	k=(i+2)%3;

	/*LST-3/9R (or CST-3/6C if alphab=0.0)*/
	*(*(Bmb+0)+3*i+0)=*(Lx+i);
	*(*(Bmb+0)+3*i+1)=0.0;
	*(*(Bmb+0)+3*i+2)=*(Lx+i)*(*(Lx+k)-*(Lx+j))*2.0*area*alphab/6.0;

	*(*(Bmb+1)+3*i+0)=0.0;
	*(*(Bmb+1)+3*i+1)=*(Ly+i);
	*(*(Bmb+1)+3*i+2)=*(Ly+i)*(*(Ly+k)-*(Ly+j))*2.0*area*alphab/6.0;

	*(*(Bmb+2)+3*i+0)=*(Ly+i);
	*(*(Bmb+2)+3*i+1)=*(Lx+i);
	*(*(Bmb+2)+3*i+2)=(*(Lx+j)**(Ly+j)-*(Lx+k)**(Ly+k))*2.0*area*alphab/3.0;
  }

  /*SHAPE FUNCTION FOR BENDING*/
  for(ii=0;ii<ngp;ii++)
  {
	/*SHAPE FUNCTION FOR IN-PLANE*/
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Q+i)+j)=1.5*sqrt(beta[0])*(*(*(L+ii)+0)**(*(Q1+i)+j)+*(*(L+ii)+1)**(*(Q2+i)+j)+*(*(L+ii)+2)**(*(Q3+i)+j));
		/*1.5 IS FOR ISOTROPIC MATERIAL*/
	  }
	}
	TnQ = matrixmatrixIII(Tn,Q,3,3,3);
	Bmh = matrixmatrixIII(TnQ,Th,3,3,9);

	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(*(B+ii)+i)+6*j+0)=*(*(Bmb+i)+3*j+0)+*(*(Bmh+i)+3*j+0);
		*(*(*(B+ii)+i)+6*j+1)=*(*(Bmb+i)+3*j+1)+*(*(Bmh+i)+3*j+1);
		*(*(*(B+ii)+i)+6*j+5)=*(*(Bmb+i)+3*j+2)+*(*(Bmh+i)+3*j+2);
	  }
	}

	freematrix(TnQ,3);
	freematrix(Bmh,3);

	/*BCIZ BY Zienkiewicz*/
	/*
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
	*/

	/*DKT:Discrete Kirchhoff Triangular BY Stricklin & Dhatt*/
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

	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(*(B+ii)+3+i)+6*j+2)=*(*(Bb+i)+3*j+0);
		*(*(*(B+ii)+3+i)+6*j+3)=*(*(Bb+i)+3*j+1);
		*(*(*(B+ii)+3+i)+6*j+4)=*(*(Bb+i)+3*j+2);
	  }
	}
  }



  freematrix(exy,3);
  freematrix(L,shell.ngp);
  free(a);
  free(b);
  free(c);
  free(Lx);
  free(Ly);
  free(len);
  free(aa);
  free(bb);
  free(cc);
  free(dd);
  free(ee);
  freematrix(Q1,3);
  freematrix(Q2,3);
  freematrix(Q3,3);
  freematrix(Th,3);
  freematrix(Tn,3);
  freematrix(Bmb,3);
  freematrix(Bb,3);

  return B;
}



double **assemshellemtx(struct oshell shell)
/*ASSEMBLAGE ELASTIC MATRIX.*/
{
  int i,j;
  int ii;/*LOOP FOR INTEGRATION POINTS*/
  int nnod = shell.nnod;
  int ngp = shell.ngp;
  int nstress = shell.nstress;
  double** Ke;
  double area;
  double* w ;
  double** C, ***B;
  double** Dm,** Db;
  double** Bm,** Bb;
  double** Km,** Kb;



  Ke=(double **)malloc(6*nnod*sizeof(double *));
  for(i=0;i<6*nnod;i++)
  {
	*(Ke+i)=(double *)malloc(6*nnod*sizeof(double));
	for(j=0;j<6*nnod;j++)
	{
	  *(*(Ke+i)+j)=0.0;                                              /*INITIAL.*/
	}
  }

  area = shell.area;
  w=(double *)malloc(ngp*sizeof(double));
  for(ii=0;ii<ngp;ii++)
  {
	*(w+ii)=area*shell.gp[ii].w;
  }

  Dm=(double **)malloc(3*sizeof(double *));
  Db=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Dm+i)=(double *)malloc(3*sizeof(double));
	*(Db+i)=(double *)malloc(3*sizeof(double));
  }
  Bm=(double **)malloc(3*sizeof(double *));
  Bb=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bm+i)=(double *)malloc(9*sizeof(double));
	*(Bb+i)=(double *)malloc(9*sizeof(double));
  }

  C = shellC(shell);
  B = shellB(shell);

  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Dm+i)+j) = *(*(C+i)+j);
	  *(*(Db+i)+j) = *(*(C+3+i)+3+j);
	}
  }

  for(ii=0;ii<ngp;ii++)
  {
	//IN-PLANE
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bm+i)+3*j+0) = *(*(*(B+ii)+i)+6*j+0);
		*(*(Bm+i)+3*j+1) = *(*(*(B+ii)+i)+6*j+1);
		*(*(Bm+i)+3*j+2) = *(*(*(B+ii)+i)+6*j+5);
	  }
	}
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

	Km = transformationEx(Dm,Bm,3,9);
	Kb = transformationEx(Db,Bb,3,9);

	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Ke+6*i+0)+6*j+0) += *(*(Km+3*i+0)+3*j+0)**(w+ii);
		*(*(Ke+6*i+0)+6*j+1) += *(*(Km+3*i+0)+3*j+1)**(w+ii);
		*(*(Ke+6*i+0)+6*j+5) += *(*(Km+3*i+0)+3*j+2)**(w+ii);
		*(*(Ke+6*i+1)+6*j+0) += *(*(Km+3*i+1)+3*j+0)**(w+ii);
		*(*(Ke+6*i+1)+6*j+1) += *(*(Km+3*i+1)+3*j+1)**(w+ii);
		*(*(Ke+6*i+1)+6*j+5) += *(*(Km+3*i+1)+3*j+2)**(w+ii);
		*(*(Ke+6*i+5)+6*j+0) += *(*(Km+3*i+2)+3*j+0)**(w+ii);
		*(*(Ke+6*i+5)+6*j+1) += *(*(Km+3*i+2)+3*j+1)**(w+ii);
		*(*(Ke+6*i+5)+6*j+5) += *(*(Km+3*i+2)+3*j+2)**(w+ii);
	  }
	}
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Ke+6*i+2)+6*j+2) += *(*(Kb+3*i+0)+3*j+0)**(w+ii);
		*(*(Ke+6*i+2)+6*j+3) += *(*(Kb+3*i+0)+3*j+1)**(w+ii);
		*(*(Ke+6*i+2)+6*j+4) += *(*(Kb+3*i+0)+3*j+2)**(w+ii);
		*(*(Ke+6*i+3)+6*j+2) += *(*(Kb+3*i+1)+3*j+0)**(w+ii);
		*(*(Ke+6*i+3)+6*j+3) += *(*(Kb+3*i+1)+3*j+1)**(w+ii);
		*(*(Ke+6*i+3)+6*j+4) += *(*(Kb+3*i+1)+3*j+2)**(w+ii);
		*(*(Ke+6*i+4)+6*j+2) += *(*(Kb+3*i+2)+3*j+0)**(w+ii);
		*(*(Ke+6*i+4)+6*j+3) += *(*(Kb+3*i+2)+3*j+1)**(w+ii);
		*(*(Ke+6*i+4)+6*j+4) += *(*(Kb+3*i+2)+3*j+2)**(w+ii);
	  }
	}

	freematrix(Km,9);
	freematrix(Kb,9);
	/*
	K = transformationEx(C,*(B+ii),nstress,6*nnod);

	for(i=0;i<6*nnod;i++)
	{
	  for(j=0;j<6*nnod;j++)
	  {
		*(*(Ke+i)+j) += *(*(K+i)+j)*(w+ii);
	  }
	}
	freematrix(K,6*nnod);
	*/
  }

  freematrix(C, nstress);
  for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
  free(B);

  freematrix(Dm,3);
  freematrix(Db,3);
  freematrix(Bm,3);
  freematrix(Bb,3);

  free(w);
  return Ke;
}


double **assemshellpmtx(struct oshell shell,double*** C,double*** B)
/*ASSEMBLAGE ELASTIC MATRIX.*/
{
  int i,j;
  int ii;/*LOOP FOR INTEGRATION POINTS*/
  int nnod = shell.nnod;
  int ngp = shell.ngp;
  int nstress = shell.nstress;
  double** Kp;
  double area;
  double* w;
  double** Dm,** Dmb,** Dbm,** Db;
  double** Bm,** Bb;
  double** Bmt,** Bbt,** DmbBb,** DbmBm;
  double** Km,** Kmb,** Kbm,** Kb;


  Kp=(double **)malloc(6*nnod*sizeof(double *));
  for(i=0;i<6*nnod;i++)
  {
	*(Kp+i)=(double *)malloc(6*nnod*sizeof(double));
	for(j=0;j<6*nnod;j++)
	{
	  *(*(Kp+i)+j)=0.0;                                              /*INITIAL.*/
	}
  }

  area = shell.area;
  w=(double *)malloc(ngp*sizeof(double));
  for(ii=0;ii<ngp;ii++)
  {
	*(w+ii)=area*shell.gp[ii].w;
  }

  Dm =(double **)malloc(3*sizeof(double *));
  Dmb=(double **)malloc(3*sizeof(double *));
  Dbm=(double **)malloc(3*sizeof(double *));
  Db =(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Dm +i)=(double *)malloc(3*sizeof(double));
	*(Dmb+i)=(double *)malloc(3*sizeof(double));
	*(Dbm+i)=(double *)malloc(3*sizeof(double));
	*(Db +i)=(double *)malloc(3*sizeof(double));
  }

  Bm=(double **)malloc(3*sizeof(double *));
  Bb=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bm+i)=(double *)malloc(9*sizeof(double));
	*(Bb+i)=(double *)malloc(9*sizeof(double));
  }

  for(ii=0;ii<ngp;ii++)
  {

	/*STRESS-STRAIN MATRIX*/
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Dm +i)+j) = *(*(*(C+ii)  +i)  +j);
		*(*(Dmb+i)+j) = *(*(*(C+ii)  +i)+3+j);
		*(*(Dbm+i)+j) = *(*(*(C+ii)+3+i)  +j);
		*(*(Db +i)+j) = *(*(*(C+ii)+3+i)+3+j);
	  }
	}
	//IN-PLANE
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bm+i)+3*j+0) = *(*(*(B+ii)+i)+6*j+0);
		*(*(Bm+i)+3*j+1) = *(*(*(B+ii)+i)+6*j+1);
		*(*(Bm+i)+3*j+2) = *(*(*(B+ii)+i)+6*j+5);
	  }
	}
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

	Km  = transformationEx(Dm,Bm,3,9);

	Bmt=matrixtransposeIII(Bm,3,9);
	DmbBb=matrixmatrixIII(Dmb,Bb,3,3,9);
	Kmb=matrixmatrixIII(Bmt,DmbBb,9,3,9);


	Bbt=matrixtransposeIII(Bb,3,9);
	DbmBm=matrixmatrixIII(Dbm,Bm,3,3,9);
	Kbm=matrixmatrixIII(Bbt,DbmBm,9,3,9);

	Kb  = transformationEx(Db,Bb,3,9);


	/*Km*/
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Kp+6*i+0)+6*j+0) += *(*(Km+3*i+0)+3*j+0)**(w+ii);
		*(*(Kp+6*i+0)+6*j+1) += *(*(Km+3*i+0)+3*j+1)**(w+ii);
		*(*(Kp+6*i+0)+6*j+5) += *(*(Km+3*i+0)+3*j+2)**(w+ii);
		*(*(Kp+6*i+1)+6*j+0) += *(*(Km+3*i+1)+3*j+0)**(w+ii);
		*(*(Kp+6*i+1)+6*j+1) += *(*(Km+3*i+1)+3*j+1)**(w+ii);
		*(*(Kp+6*i+1)+6*j+5) += *(*(Km+3*i+1)+3*j+2)**(w+ii);
		*(*(Kp+6*i+5)+6*j+0) += *(*(Km+3*i+2)+3*j+0)**(w+ii);
		*(*(Kp+6*i+5)+6*j+1) += *(*(Km+3*i+2)+3*j+1)**(w+ii);
		*(*(Kp+6*i+5)+6*j+5) += *(*(Km+3*i+2)+3*j+2)**(w+ii);
	  }
	}
	/*Kmb*/
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Kp+6*i+0)+6*j+2) += *(*(Kmb+3*i+0)+3*j+0)**(w+ii);
		*(*(Kp+6*i+0)+6*j+3) += *(*(Kmb+3*i+0)+3*j+1)**(w+ii);
		*(*(Kp+6*i+0)+6*j+4) += *(*(Kmb+3*i+0)+3*j+2)**(w+ii);
		*(*(Kp+6*i+1)+6*j+2) += *(*(Kmb+3*i+1)+3*j+0)**(w+ii);
		*(*(Kp+6*i+1)+6*j+3) += *(*(Kmb+3*i+1)+3*j+1)**(w+ii);
		*(*(Kp+6*i+1)+6*j+4) += *(*(Kmb+3*i+1)+3*j+2)**(w+ii);
		*(*(Kp+6*i+5)+6*j+2) += *(*(Kmb+3*i+2)+3*j+0)**(w+ii);
		*(*(Kp+6*i+5)+6*j+3) += *(*(Kmb+3*i+2)+3*j+1)**(w+ii);
		*(*(Kp+6*i+5)+6*j+4) += *(*(Kmb+3*i+2)+3*j+2)**(w+ii);
	  }
	}
	/*Kbm*/
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Kp+6*i+2)+6*j+0) += *(*(Kbm+3*i+0)+3*j+0)**(w+ii);
		*(*(Kp+6*i+2)+6*j+1) += *(*(Kbm+3*i+0)+3*j+1)**(w+ii);
		*(*(Kp+6*i+2)+6*j+5) += *(*(Kbm+3*i+0)+3*j+2)**(w+ii);
		*(*(Kp+6*i+3)+6*j+0) += *(*(Kbm+3*i+1)+3*j+0)**(w+ii);
		*(*(Kp+6*i+3)+6*j+1) += *(*(Kbm+3*i+1)+3*j+1)**(w+ii);
		*(*(Kp+6*i+3)+6*j+5) += *(*(Kbm+3*i+1)+3*j+2)**(w+ii);
		*(*(Kp+6*i+4)+6*j+0) += *(*(Kbm+3*i+2)+3*j+0)**(w+ii);
		*(*(Kp+6*i+4)+6*j+1) += *(*(Kbm+3*i+2)+3*j+1)**(w+ii);
		*(*(Kp+6*i+4)+6*j+5) += *(*(Kbm+3*i+2)+3*j+2)**(w+ii);
	  }
	}
	/*Kb*/
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Kp+6*i+2)+6*j+2) += *(*(Kb+3*i+0)+3*j+0)**(w+ii);
		*(*(Kp+6*i+2)+6*j+3) += *(*(Kb+3*i+0)+3*j+1)**(w+ii);
		*(*(Kp+6*i+2)+6*j+4) += *(*(Kb+3*i+0)+3*j+2)**(w+ii);
		*(*(Kp+6*i+3)+6*j+2) += *(*(Kb+3*i+1)+3*j+0)**(w+ii);
		*(*(Kp+6*i+3)+6*j+3) += *(*(Kb+3*i+1)+3*j+1)**(w+ii);
		*(*(Kp+6*i+3)+6*j+4) += *(*(Kb+3*i+1)+3*j+2)**(w+ii);
		*(*(Kp+6*i+4)+6*j+2) += *(*(Kb+3*i+2)+3*j+0)**(w+ii);
		*(*(Kp+6*i+4)+6*j+3) += *(*(Kb+3*i+2)+3*j+1)**(w+ii);
		*(*(Kp+6*i+4)+6*j+4) += *(*(Kb+3*i+2)+3*j+2)**(w+ii);
	  }
	}
	freematrix(Km,3*nnod);
	freematrix(Kmb,3*nnod);
	freematrix(Kbm,3*nnod);
	freematrix(Kb,3*nnod);
	freematrix(Bmt,9);
	freematrix(Bbt,9);
	freematrix(DmbBb,3);
	freematrix(DbmBm,3);


	/*
	K = transformationEx(*(C+ii),*(B+ii),nstress,6*nnod);

	for(i=0;i<6*nnod;i++)
	{
	  for(j=0;j<6*nnod;j++)
	  {
		*(*(Kp+i)+j) += *(*(K+i)+j)*(w+ii);
	  }
	}
	freematrix(K,6*nnod);
	*/
  }

  freematrix(Dm,3);
  freematrix(Db,3);
  freematrix(Bm,3);
  freematrix(Bb,3);

  free(w);
  return Kp;
}




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



	//DRILL
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

double **assemshellmmtx(struct oshell shell)
/*ASSEMBLAGE ELASTIC MATRIX.*/
{
  char string[100];
  int nnod = shell.nnod;
  int ngp = shell.ngp;
  int i,j,k,ii;
  double **m,**exy,**N,**Nt,**NtN,**L;
  double t,hiju;
  double *b,*c;
  double *w;
  double area;
  double Liij,Lijj,Liik,Likk,Lijk;
  double alpha=0.0;
  //double alpha=0.15;

  /*

	mvct = (double *)malloc(msize * sizeof(double));
	for (i = 0; i < msize; i++)
	{
		*(mvct + i) = 0.0;
	}

	for (i = 1; i <= nshell; i++)
	{
		Me = transformationIII(Me, T, 6*nnod);

		masstotal = 0.0;
		massdiag = 0.0;
		for (ii = 0; ii < 18; ii++)
		{
			for (jj = 0; jj < 18; jj++)
			{
				masstotal += Me[ii][jj];
				if (ii == jj && ii % 6 < 3)massdiag += Me[ii][jj];
			}
		}
		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				*(mvct + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(*(Me + 6 * ii + jj) + 6 * ii + jj) * masstotal / massdiag;
			}
		}
	}

  */

  t=shell.sect->area;
  hiju=shell.sect->hiju[0];

  m=(double **)malloc(6*nnod*sizeof(double *));
  for(i=0;i<6*nnod;i++)
  {
	*(m+i)=(double *)malloc(6*nnod*sizeof(double));
	for(j=0;j<6*nnod;j++)
	{
	  *(*(m+i)+j)=0.0;                                              /*INITIAL.*/
	}
  }
  area = shell.area;


#if 0
  /*WEIGHT OF EACH NODE*/

  w=(double *)malloc(ngp*sizeof(double));
  for(ii=0;ii<ngp;ii++)
  {
	*(w+ii)=area*shell.gp[ii].w;
  }

  exy=shelllocalcoord(shell);/*shell local coordinate {xi,yi(,zi=0)}*/
  /*TRIANGLE COORDINATION*/
  L=(double **)malloc(ngp*sizeof(double *));
  for(i=0;i<ngp;i++)
  {
	*(L+i)=(double *)malloc(3*sizeof(double));
  }
  /*AREA COORD Li=(ai+bix+ciy)/(2*area)*/
  if(nnod==3 && ngp==7)
  {
	  /*3 NODES 7 GAUSS POINTS*/
	  *(*(L+0)+0)=1.0/3.0;*(*(L+0)+1)=1.0/3.0;*(*(L+0)+2)=1.0/3.0;

	  *(*(L+1)+0)=1.0/2.0;*(*(L+1)+1)=1.0/2.0;*(*(L+1)+2)=0.0;
	  *(*(L+2)+0)=0.0;*(*(L+2)+1)=1.0/2.0;*(*(L+2)+2)=1.0/2.0;
	  *(*(L+3)+0)=1.0/2.0;*(*(L+3)+1)=0.0;*(*(L+3)+2)=1.0/2.0;

	  *(*(L+4)+0)=1.0;*(*(L+4)+1)=0.0;*(*(L+4)+2)=0.0;
	  *(*(L+5)+0)=0.0;*(*(L+5)+1)=1.0;*(*(L+5)+2)=0.0;
	  *(*(L+6)+0)=0.0;*(*(L+6)+1)=0.0;*(*(L+6)+2)=1.0;
  }
  if(nnod==3 && ngp==3)
  {
	  /*3 NODES 3 GAUSS POINTS*/
	  *(*(L+0)+0)=1.0/6.0;*(*(L+0)+1)=1.0/6.0;*(*(L+0)+2)=2.0/3.0;
	  *(*(L+1)+0)=1.0/6.0;*(*(L+1)+1)=2.0/3.0;*(*(L+1)+2)=1.0/6.0;
	  *(*(L+2)+0)=2.0/3.0;*(*(L+2)+1)=1.0/6.0;*(*(L+2)+2)=1.0/6.0;
  }
  if(nnod==3 && ngp==4)
  {
	  /*3 NODES 4 GAUSS POINTS*/
	  *(*(L+0)+0)=1.0/3.0;*(*(L+0)+1)=1.0/3.0;*(*(L+0)+2)=1.0/3.0;
	  *(*(L+1)+0)=1.0/6.0;*(*(L+1)+1)=1.0/6.0;*(*(L+1)+2)=2.0/3.0;
	  *(*(L+2)+0)=1.0/6.0;*(*(L+2)+1)=2.0/3.0;*(*(L+2)+2)=1.0/6.0;
	  *(*(L+3)+0)=2.0/3.0;*(*(L+3)+1)=1.0/6.0;*(*(L+3)+2)=1.0/6.0;
  }

  b=(double *)malloc(3*sizeof(double));
  c=(double *)malloc(3*sizeof(double));
  for(i=0;i<nnod;i++)
  {
	  j=i+1;
	  k=i+2;
	  if(j>=3)j-=3;
	  if(k>=3)k-=3;
	  *(b+i)=*(*(exy+j)+1)-*(*(exy+k)+1);/*yj-yk*/
	  *(c+i)=*(*(exy+k)+0)-*(*(exy+j)+0);/*xk-xj*/
  }



  for(ii=0;ii<ngp;ii++)
  {
	N=(double **)malloc(nnod*sizeof(double *));
	for(i=0;i<nnod;i++)
	{
	  *(N+i)=(double *)malloc(6*nnod*sizeof(double));
	  for(j=0;j<6*nnod;j++)
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
	Nt=matrixtransposeIII(N,nnod,6*nnod);
	NtN=matrixmatrixIII(Nt,N,6*nnod,nnod,6*nnod);
	for(i=0;i<6*nnod;i++)
	{
	  for(j=0;j<6*nnod;j++)
	  {
		*(*(m+i)+j) += *(*(NtN+i)+j)**(w+ii)*t*hiju;
	  }
	}
	freematrix(N,nnod);
	freematrix(Nt,6*nnod);
	freematrix(NtN,6*nnod);
  }

  free(w);
  freematrix(L,ngp);
  freematrix(exy,3);
  free(b);
  free(c);
#endif
#if 1
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  for(k=0;k<3;k++)
	  {
		if(i==j)
		{
		  *(*(m+6*i+k+0)+6*j+k+0) = hiju*area*t/6.0;
		  *(*(m+6*i+k+3)+6*j+k+3) = hiju*area*pow(t,3)/18.0;
		}
		else
		{
		  *(*(m+6*i+k+0)+6*j+k+0) = hiju*area*t/12.0;
		  *(*(m+6*i+k+3)+6*j+k+3) = hiju*area*pow(t,3)/36.0;
		}
	  }
	}
  }
  /*MASS MATRIX BY ZHONG AND ALMEDIA.*/
#endif



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
  /*AREA COORD Li=(ai+bix+ciy)/(2*area)*/
  //if(nnod==3 && ngp==7)
  //{
	  /*3 NODES 7 GAUSS POINTS*/
	  *(*(L+0)+0)=1.0/3.0;*(*(L+0)+1)=1.0/3.0;*(*(L+0)+2)=1.0/3.0;

	  *(*(L+1)+0)=1.0/2.0;*(*(L+1)+1)=1.0/2.0;*(*(L+1)+2)=0.0;
	  *(*(L+2)+0)=0.0;*(*(L+2)+1)=1.0/2.0;*(*(L+2)+2)=1.0/2.0;
	  *(*(L+3)+0)=1.0/2.0;*(*(L+3)+1)=0.0;*(*(L+3)+2)=1.0/2.0;

	  *(*(L+4)+0)=1.0;*(*(L+4)+1)=0.0;*(*(L+4)+2)=0.0;
	  *(*(L+5)+0)=0.0;*(*(L+5)+1)=1.0;*(*(L+5)+2)=0.0;
	  *(*(L+6)+0)=0.0;*(*(L+6)+1)=0.0;*(*(L+6)+2)=1.0;
  //}
  /*
  if(nnod==3 && ngp==3)
  {
	  *(*(L+0)+0)=1.0/6.0;*(*(L+0)+1)=1.0/6.0;*(*(L+0)+2)=2.0/3.0;
	  *(*(L+1)+0)=1.0/6.0;*(*(L+1)+1)=2.0/3.0;*(*(L+1)+2)=1.0/6.0;
	  *(*(L+2)+0)=2.0/3.0;*(*(L+2)+1)=1.0/6.0;*(*(L+2)+2)=1.0/6.0;
  }
  if(nnod==3 && ngp==4)
  {
	  *(*(L+0)+0)=1.0/3.0;*(*(L+0)+1)=1.0/3.0;*(*(L+0)+2)=1.0/3.0;
	  *(*(L+1)+0)=1.0/6.0;*(*(L+1)+1)=1.0/6.0;*(*(L+1)+2)=2.0/3.0;
	  *(*(L+2)+0)=1.0/6.0;*(*(L+2)+1)=2.0/3.0;*(*(L+2)+2)=1.0/6.0;
	  *(*(L+3)+0)=2.0/3.0;*(*(L+3)+1)=1.0/6.0;*(*(L+3)+2)=1.0/6.0;
  }
  */


  b=(double *)malloc(3*sizeof(double));
  c=(double *)malloc(3*sizeof(double));
  q=(double *)malloc(3*sizeof(double));
  exy=shelllocalcoord(shell);
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

		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell);

		*volume += shellvolume(shell, drccos);                   		/*VOLUME*/

		freematrix(drccos, 3);
	}
	return;
}

/*shells or mshell->shell*/
void inputshell(struct oshell *shells,
				struct memoryshell *mshell,int offset,
				struct oshell *shell)
{
  int i,j,k;
  shell->nnod=(shells+offset)->nnod;
  shell->ngp=(shells+offset)->ngp;
  shell->nstress=(shells+offset)->nstress;
  shell->loff=offset;
  shell->code=(shells+offset)->code;                  /*ELEMENT CODE.*/
  shell->sect=(shells+offset)->sect;               /*SECTION POINTER.*/
  shell->prate=(shells+offset)->prate;
  shell->brate=(shells+offset)->brate;
  shell->area=(shells+offset)->area;



  if(mshell!=NULL)
  {
	  for(i=0;i<shell->nnod;i++)
	  {
		  shell->node[i]=(shells+offset)->node[i];
		  for(j=0;j<6;j++)                                      /*STRESS.*/
		  {
			shell->stress[i][j]=(mshell+offset)->stress[i][j];
		  }
	  }
	  memcpy(shell->gp,(mshell+offset)->gp,sizeof(shell->gp));
  }
  else
  {
     for(i=0;i<shell->nnod;i++)
	  {
		  shell->node[i]=(shells+offset)->node[i];
		  for(j=0;j<6;j++)                                      /*STRESS.*/
		  {
			shell->stress[i][j]=(shells+offset)->stress[i][j];
		  }
	  }
	  memcpy(shell->gp,(shells+offset)->gp,sizeof(shell->gp));
  }
  return;
}/*inputshell*/

/*shell -> shells*/
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
  memcpy((shells+offset)->gp,shell->gp,sizeof(shell->gp));
  return;
}/*inputshell*/

/*shells -> mshell*/
void initialshell(struct oshell *shells,struct memoryshell *mshell,int nshell)
{
  int i,j,k;

  for(i = 0; i < nshell; i++)
  {
	  (mshell+i)->code=(shells+i)->code;

	  for(j=0;j<(shells+i)->nnod;j++)
	  {
		for(k=0;k<6;k++)                                      /*STRESS.*/
		{
		  (mshell+i)->stress[j][k]=(shells+i)->stress[j][k];
		}
	  }
	  memcpy((mshell+i)->gp,(shells+i)->gp,sizeof((shells+i)->gp));

	  for(j=0;j<(shells+i)->ngp;j++)
	  {
		((mshell+i)->gp[j]).lambda[0]=0.0/*((shells+i)->gp[j]).lambda[0]*/;
		((mshell+i)->gp[j]).lambda[1]=0.0/*((shells+i)->gp[j]).lambda[1]*/;

		for(k=0;k<(shells+i)->gp[j].ngp2;k++)
		{
			((mshell+i)->gp[j]).gp2[j].lambda=0.0/*((shells+i)->gp[j]).gp2[j].lambda*/;
		}
	  }
  }
  return;
}/*initialshell*/


/*mshell->shell*/




double yieldstress(struct osect* sect, double alpha, double* dy, double* ddy)
{
	char str[256];
	double y;
	double yinit,c,n;
	yinit = sect->yieldinit;
	c = sect->yieldcoefficient;
	n = sect->yieldpower;

	double sa,sc,sn;/*SWIFT'S PARAMS*/

	if(n==0.0)
	{
	  y = yinit;
	  if(dy!=NULL)*dy = 0.0;
	  if(ddy!=NULL)*ddy = 0.0;
	}
	/*
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
	*/
	else
	{
	  /*SWIFT'S HARDNING LAW : σ=sc*(sa+ε_p)^sn*/
	  /*RAMBERG-OSGOOD : ε=σ/E+(σ/c)^n, ε_p=(σ/c)^n*/
	  /*MODIFIED RAMBERG-OSGOOD : ε=σ/E+(σ/c)^n+(σy/c)^n, ε_p=(σ/c)^n+(σy/c)^n, σ=c{(σy/c)^n+ε_p}^(1/n)*/

	  sa = pow(yinit/c,n);/*(σy/c)^n*/
	  sc = c;
	  sn = 1/n;

	  y = sc * pow(sa+alpha,sn);
	  if(sn==1.0)
	  {
		if(dy!=NULL)*dy = sc;
		if(ddy!=NULL)*ddy = 0.0;
	  }
	  else
	  {
		if(dy!=NULL)*dy = sc*sn*pow(sa+alpha,n-1);
		if(ddy!=NULL)*ddy = sc*sn*sn*pow(sa+alpha,n-2);
	  }
	}




	return y;
}


double*** shellCconsistent(struct oshell shell)
{
  int ii;
  double *** consistentC;/*ngp x nstress x nstress*/

  consistentC = (double***)malloc(shell.ngp * sizeof(double**));

  for(ii = 0; ii < shell.ngp; ii++)/*FOR EACH INTEGRATION POINT*/
  {
	if(shell.gp[ii].ngp2==0)
	{
		*(consistentC+ii) = shellCconsistentilyushin(shell, ii);

	}
	else
	{
		*(consistentC+ii) = shellCconsistentmises(shell, ii);
    }
  }
  return consistentC;
}


/*ELASTO-PLASTICITY*/
double* ilyushin(struct oshell* shell, int ii, double* lambda, int UPDATEFLAG)
{
  int i,j;
  int nstress = shell->nstress;
  struct gausspoint* gp;
  double* qstress,* ostress;
  double* f;
  double fn,fm;
  double fn2,fm2,fnfm;
  double E,t,poi;
  double yinit,y;
  double alpha;
  double qn,qm,qnm;

  double v[3],p[3],e[3],det[3];
  double O[2][2][3], Oinv[2][2][3];

  double c[3],a[6];
  double ** g;

  char str[800];

  gp = &(shell->gp[ii]);

  qstress=(double*)malloc(nstress * sizeof(double));
  ostress=(double*)malloc(nstress * sizeof(double));
  f=(double*)malloc(2 * sizeof(double));

  /*          | 1 -1  0 |  */
  /*[Q] = 1/√2| 1  1  0 |  */
  /*          | 0  0 √2 |  */

  /*[Q]diag[a,b,c][Q]^T  */
  /*    | a+b a-b   0 |  */
  /*=1/2| a-b a+b   0 |  */
  /*    |   0   0  2c |  */

  /*{σtry-αtry}=([I]+2Δλ[C][A]){σ-α}=[W]^-1{σ-α} */
  /*{σ-α}=[W]{σtry-αtry}                         */
  /*[W]^-1=(diag[Q,Q])[O](diag[Q^T,Q^T])         */
  /*[W]=(diag[Q,Q])[O]^-1(diag[Q^T,Q^T])         */

  E=shell->sect->E;
  t=shell->sect->area;
  poi=shell->sect->poi;
  yinit=shell->sect->yieldinit;

  fn=shell->sect->fmax[0] /*yinit*t*/;
  fm=shell->sect->fmax[3] /*yinit*pow(t,2)/4.0*/;

  fn2=pow(fn,2);
  fm2=pow(fm,2);
  fnfm=fn*fm;


  v[0]=E/(1-poi);
  v[1]=E/(1+poi);
  v[2]=E/(2.0*(1+poi));

  p[0]=0.5;
  p[1]=1.5;
  p[2]=3.0;

  e[0] = (*(lambda+0) + *(lambda+1)) / fn2;
  e[1] = (*(lambda+0) + *(lambda+1)) / fm2;
  e[2] = (*(lambda+0) - *(lambda+1)) /(2.0*sqrt(3.0)*fnfm);


  /*[O]*/
  for (i = 0; i < 3; i++)
  {
	O[0][0][i] = 2.0*t            *v[i]*p[i] /*+ 4.0/3.0*Hkin*[p]*/; O[0][1][i] = O[0][0][i];
	O[1][0][i] = 2.0*pow(t,3)/12.0*v[i]*p[i] /*+ 4.0/3.0*Hkin*[p]*/; O[1][1][i] = O[1][0][i];
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




  /*(diag[Q^T,Q^T]){σtry-αtry}*/
  *(qstress+0)=( gp->stress[0] + gp->stress[1] )/sqrt(2.0);
  *(qstress+1)=(-gp->stress[0] + gp->stress[1] )/sqrt(2.0);
  *(qstress+2)=  gp->stress[2];
  *(qstress+3)=( gp->stress[3] + gp->stress[4] )/sqrt(2.0);
  *(qstress+4)=(-gp->stress[3] + gp->stress[4] )/sqrt(2.0);
  *(qstress+5)=  gp->stress[5];

  /*[O]^-1*(diag[Q^T,Q^T]){σtry-αtry}*/
  *(ostress+0) = Oinv[0][0][0]**(qstress+0) + Oinv[0][1][0]**(qstress+3);
  *(ostress+1) = Oinv[0][0][1]**(qstress+1) + Oinv[0][1][1]**(qstress+4);
  *(ostress+2) = Oinv[0][0][2]**(qstress+2) + Oinv[0][1][2]**(qstress+5);
  *(ostress+3) = Oinv[1][0][0]**(qstress+0) + Oinv[1][1][0]**(qstress+3);
  *(ostress+4) = Oinv[1][0][1]**(qstress+1) + Oinv[1][1][1]**(qstress+4);
  *(ostress+5) = Oinv[1][0][2]**(qstress+2) + Oinv[1][1][2]**(qstress+5);



  /*f={σtry-αtry}^T[W]^T[A][W]{σtry-αtry}*/

  qn = (p[0]**(ostress+0)**(ostress+0) + p[1]**(ostress+1)**(ostress+1) + p[2]**(ostress+2)**(ostress+2))/fn2 ;
  qm = (p[0]**(ostress+3)**(ostress+3) + p[1]**(ostress+4)**(ostress+4) + p[2]**(ostress+5)**(ostress+5))/fm2 ;
  qnm= (p[0]**(ostress+0)**(ostress+3) + p[1]**(ostress+1)**(ostress+4) + p[2]**(ostress+2)**(ostress+5))/fnfm;

  *(f+0) = qn+qm+(qnm/sqrt(3.0));
  *(f+1) = qn+qm-(qnm/sqrt(3.0));

  /*
  qn  = ( pow(gp->stress[0],2)+pow(gp->stress[1],2)-(gp->stress[0])*(gp->stress[1])+3.0*pow(gp->stress[2],2) )/fn2;
  qm  = ( pow(gp->stress[3],2)+pow(gp->stress[4],2)-(gp->stress[3])*(gp->stress[4])+3.0*pow(gp->stress[5],2) )/fm2;
  qnm = (       (gp->stress[0])*(gp->stress[3])
		  +     (gp->stress[1])*(gp->stress[4])
		  -0.5* (gp->stress[0])*(gp->stress[4])
		  -0.5* (gp->stress[1])*(gp->stress[3])
		  +3.0* (gp->stress[2])*(gp->stress[5])  )/fnfm;




  *(f+0) = p[0] * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 + *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + p[1] * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 + *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + p[2] * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 + *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm));
  *(f+1) = p[0] * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 - *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + p[1] * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 - *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + p[2] * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 - *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm));

  */


  alpha = gp->alpha + 2.0*( *(lambda+0)*sqrt(*(f+0)) + *(lambda+1)*sqrt(*(f+1)) )/(t*yinit);
  y = yieldstress(shell->sect,alpha,NULL,NULL);
  *(f+0) -= pow(y/yinit,2.0);
  *(f+1) -= pow(y/yinit,2.0);

  if(UPDATEFLAG == 1)/*UPDATE PARAMS*/
  {
	  /*(diag[Q,Q])[O]^-1*(diag[Q^T,Q^T]){σtry-αtry}*/
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

	  gp->f[0] = *(f+0);
	  gp->f[1] = *(f+1);

	  gp->lambda[0] = *(lambda+0);
	  gp->lambda[1] = *(lambda+1);

	  c[0] = fn*fn;
	  c[1] = fm*fm;
	  c[2] = 2.0*sqrt(3.0)*fn*fm;

	  a[0] = 2.0 * gp->stress[0] - gp->stress[1];
	  a[1] = 2.0 * gp->stress[1] - gp->stress[0];
	  a[2] = 6.0 * gp->stress[2];
	  a[3] = 2.0 * gp->stress[3] - gp->stress[4];
	  a[4] = 2.0 * gp->stress[4] - gp->stress[3];
	  a[5] = 6.0 * gp->stress[5];

	  g = (double**)malloc(2 * sizeof(double*));
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

	  *(*(g+1)+0) =   a[0]/c[0] - a[3]/c[2];
	  *(*(g+1)+1) =   a[1]/c[0] - a[4]/c[2];
	  *(*(g+1)+2) =   a[2]/c[0] - a[5]/c[2];
	  *(*(g+1)+3) = - a[0]/c[2] + a[3]/c[1];
	  *(*(g+1)+4) = - a[1]/c[2] + a[4]/c[1];
	  *(*(g+1)+5) = - a[2]/c[2] + a[5]/c[1];

	  /*UPDATE HARDENING PARAMS*/
	  for (i = 0; i < nstress; i++)
	  {
		  gp->estrain[i] -= *(lambda+0)**(*(g+0)+i) + *(lambda+1)**(*(g+1)+i);
		  gp->pstrain[i] += *(lambda+0)**(*(g+0)+i) + *(lambda+1)**(*(g+1)+i);
	  }
	  freematrix(g,2);

  }
  free(qstress);
  free(ostress);

  return f;
}



double** shellCconsistentilyushin(struct oshell shell, int ii)
{
  int i,j;
  int nnod,ngp,nstress;

  struct gausspoint* gp;

  double fn,fm;
  double fn2,fm2,fnfm;
  double E,t,A,I,poi;
  double v[3],p[3],e[3],det[3];
  double O[2][2][3], Oinv[2][2][3];
  double c[3],a[6];

  double ** C, ** consistentC;
  double ** W;
  double ** H, ** g;
  double Halpha, galpha;
  double y, dy, ddy;
  double yinit;

  double* g0,* g1,* N0,* N1;
  double M00,M01,M10,M11;
  double detM;
  double beta;
  double tolerance = 1.0e-6;/*SAME AS RETURNMAP*/

  char str[800];



  nnod = shell.nnod;
  ngp = shell.ngp;
  nstress = shell.nstress;

  E=shell.sect->E;
  t=shell.sect->area;
  poi=shell.sect->poi;
  yinit = shell.sect->yieldinit;

  fn=shell.sect->fmax[0] /*yinit*t*/;
  fm=shell.sect->fmax[3] /*yinit*pow(t,2)/4.0*/;

  fn2=pow(fn,2);
  fm2=pow(fm,2);
  fnfm=fn*fm;

  v[0]=E/(1-poi);
  v[1]=E/(1+poi);
  v[2]=E/(2.0*(1+poi));

  p[0]=0.5;
  p[1]=1.5;
  p[2]=3.0;

  c[0] = fn*fn;
  c[1] = fm*fm;
  c[2] = 2.0*sqrt(3.0)*fn*fm;





	consistentC = (double**)malloc(nstress * sizeof(double*));
	for (i = 0; i < nstress; i++)
	{
		*(consistentC + i) = (double*)malloc(nstress * sizeof(double));
		for (j = 0; j < nstress; j++)
		{
			*(*(consistentC + i) + j) = 0.0;
		}
	}

	W = (double**)malloc(nstress * sizeof(double*));
	for (i = 0; i < nstress; i++)
	{
	  *(W + i) = (double*)malloc(nstress * sizeof(double));
	}
	g  = (double**)malloc(2 * sizeof(double*));
	for (i = 0; i < 2; i++)
	{
	  *(g+i) = (double*)malloc(nstress * sizeof(double));
	}

	gp = &(shell.gp[ii]);

	C = shellC(shell);

	e[0] = (gp->lambda[0] + gp->lambda[1]) / fn2;
	e[1] = (gp->lambda[0] + gp->lambda[1]) / fm2;
	e[2] = (gp->lambda[0] - gp->lambda[1]) /(2.0*sqrt(3.0)*fnfm);


	/*[O]*/
	for (i = 0; i < 3; i++)
	{
	  O[0][0][i] = 2.0*            t*v[i]*p[i] /*+ 4.0/3.0*Hkin*[p]*/; O[0][1][i] = O[0][0][i];
	  O[1][0][i] = 2.0*pow(t,3)/12.0*v[i]*p[i] /*+ 4.0/3.0*Hkin*[p]*/; O[1][1][i] = O[1][0][i];
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

	/*[W]=[Q][O]^-1[Q]*/
	for (i = 0; i < 2; i++)
	{
	  for (j = 0; j < 2; j++)
	  {
		  *(*(W+3*i+0)+3*j+0) = 0.5*(Oinv[i][j][0]+Oinv[i][j][1]);
		  *(*(W+3*i+0)+3*j+1) = 0.5*(Oinv[i][j][0]-Oinv[i][j][1]);
		  *(*(W+3*i+0)+3*j+2) = 0.0;
		  *(*(W+3*i+1)+3*j+0) = *(*(W+3*i+0)+3*j+1);
		  *(*(W+3*i+1)+3*j+1) = *(*(W+3*i+0)+3*j+0);
		  *(*(W+3*i+1)+3*j+2) = 0.0;
		  *(*(W+3*i+2)+3*j+0) = 0.0;
		  *(*(W+3*i+2)+3*j+1) = 0.0;
		  *(*(W+3*i+2)+3*j+2) = Oinv[i][j][2];
	  }
	}
	H = matrixmatrix(W,C,nstress);

	/*{df/dσ}={d(σAσ)/dσ}=2[A]{σ}={g}*/
	a[0] = 2.0 * gp->stress[0] - gp->stress[1];
	a[1] = 2.0 * gp->stress[1] - gp->stress[0];
	a[2] = 6.0 * gp->stress[2];
	a[3] = 2.0 * gp->stress[3] - gp->stress[4];
	a[4] = 2.0 * gp->stress[4] - gp->stress[3];
	a[5] = 6.0 * gp->stress[5];


	*(*(g+0)+0) = a[0]/c[0] + a[3]/c[2];
	*(*(g+0)+1) = a[1]/c[0] + a[4]/c[2];
	*(*(g+0)+2) = a[2]/c[0] + a[5]/c[2];
	*(*(g+0)+3) = a[0]/c[2] + a[3]/c[1];
	*(*(g+0)+4) = a[1]/c[2] + a[4]/c[1];
	*(*(g+0)+5) = a[2]/c[2] + a[5]/c[1];

	*(*(g+1)+0) =   a[0]/c[0] - a[3]/c[2];
	*(*(g+1)+1) =   a[1]/c[0] - a[4]/c[2];
	*(*(g+1)+2) =   a[2]/c[0] - a[5]/c[2];
	*(*(g+1)+3) = - a[0]/c[2] + a[3]/c[1];
	*(*(g+1)+4) = - a[1]/c[2] + a[4]/c[1];
	*(*(g+1)+5) = - a[2]/c[2] + a[5]/c[1];


	/*df/dq=d/dq(-σy^2/σyinit^2)=-2*σy*σy'(q)/(σyinit^2)=beta*/
	y = yieldstress(shell.sect, gp->alpha, &dy, &ddy);
	if(dy!=0)
	{
		Halpha = 1.0/(t*dy-((gp->lambda[0])+(gp->lambda[1]))*2.0*(dy*dy+y*ddy)/pow(yinit,2.0));
		galpha = -2.0*y*dy/pow(yinit,2.0);

		beta = Halpha*pow(galpha,2.0);
	}
	else
	{
		beta = 0;
    }

	/*CONSISTENT STIFFNESS*/
	if(abs(gp->f[0])<tolerance && abs(gp->f[1])<tolerance)
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
	  detM = M00*M11-M01*M10;

	  for (i = 0; i < nstress; i++)
	  {
		for (j = 0; j < nstress; j++)
		{
		  *(*(consistentC + i) + j) = *(*(H+i)+j) - ( M11**(N0+i)**(N0+j) - M01**(N0+i)**(N1+j) - M10**(N1+i)**(N0+j) + M00**(N1+i)**(N1+j) )/detM;
		}
	  }
	  free(N0);
	  free(N1);
	}
	else if(abs(gp->f[0])<tolerance)
	{
	  N0=matrixvector(H,*(g+0),nstress);
	  M00=dotproduct(*(g+0),N0,nstress);
	  M00+=beta;

	  for (i = 0; i < nstress; i++)
	  {
		for (j = 0; j < nstress; j++)
		{
		  *(*(consistentC + i) + j) = *(*(H+i)+j) - *(N0+i)**(N0+j)/M00;
		}
	  }
	  free(N0);
	}
	else if(abs(gp->f[1])<tolerance)
	{
	  N1=matrixvector(H,*(g+1),nstress);
	  M11=dotproduct(*(g+1),N1,nstress);
	  M11+=beta;

	  for (i = 0; i < nstress; i++)
	  {
		for (j = 0; j < nstress; j++)
		{
		  *(*(consistentC + i) + j) = *(*(H+i)+j) - *(N1+i)**(N1+j)/M11;
		}
	  }
	  free(N1);
	}
	else
	{
	  for (i = 0; i < nstress; i++)
	  {
		for (j = 0; j < nstress; j++)
		{
		  *(*(consistentC + i) + j) = *(*(H+i)+j);
		}
	  }
	}
	freematrix(H,nstress);

	freematrix(W,nstress);
	freematrix(g,2);
	freematrix(C, nstress);

	return consistentC;
}






void returnmapilyushin(struct oshell* shell, int ii)
{
   int i,j;
   //struct gausspoint* gp;
   double* lambda,* lambda01,* lambda0,* lambda1,* lambdaeps,* dlambda;
   double* f,* feps,* finit,* f01,* f0,* f1;
   double** dfdl;
   double det;
   double residual;
   double tolerance = 1.0e-6;/*SAME AS CONSISTENT_C*/
   double eps = 1.0e-2;
   char str[800];
   int iteration;
   int maxiteration=20;
   int FLAG=0;



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
   f=ilyushin(shell, ii, lambda, NULL);/*ILYUSHIN'S YIELD FUNCTION.{σ-α}^T[A]{σ-α}*/
   *(finit + 0) = *(f + 0);
   *(finit + 1) = *(f + 1);


   if(*(f+0)>0 || *(f+1)>0)/*RETURN-MAPPING PROCEDURE USING NEWTON-RAPTHON METHOD FOR ILYUSHIN'S YIELD CONDITION.*/
   {

	   /*BOTH YIELD SURFACES ARE ACTIVE.*/
	   *(lambda + 0) = 0.0;
	   *(lambda + 1) = 0.0;

	   *(f + 0) = *(finit + 0);
	   *(f + 1) = *(finit + 1);

	   residual=1.0;
	   iteration=1;
	   while(residual>tolerance && iteration<maxiteration)
	   {
		 /*SENSITIVITY*/
		 for(i=0;i<2;i++)
		 {
		   for(j=0;j<2;j++)
		   {
			 if(i==j)
			 {
				*(lambdaeps+j)=*(lambda+j)+eps;
			 }
			 else
			 {
				*(lambdaeps+j)=*(lambda+j);
			 }
		   }
		   feps=ilyushin(shell, ii, lambdaeps, NULL);


		   for(j=0;j<2;j++)
		   {
			 *(*(dfdl+j)+i)=(*(feps+j)-*(f+j))/eps;
		   }
		   free(feps);
		 }

		 /*UPDATE*/
		 det=*(*(dfdl+0)+0)* *(*(dfdl+1)+1) - *(*(dfdl+0)+1)* *(*(dfdl+1)+0);
		 if(det==0.0)
		 {
		   FLAG++;
		   break;
		 }

		 *(dlambda + 0) = -( *(*(dfdl+1)+1)**(f+0)-*(*(dfdl+0)+1)**(f+1))/det;
		 *(dlambda + 1) = -(-*(*(dfdl+1)+0)**(f+0)+*(*(dfdl+0)+0)**(f+1))/det;

		 *(lambda + 0) += *(dlambda + 0);
		 *(lambda + 1) += *(dlambda + 1);

		 /*NEW JUDGE*/
		 free(f);
		 f=ilyushin(shell, ii, lambda, NULL);

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
	   while(residual>tolerance && iteration<maxiteration)
	   {
		 /*SENSITIVITY*/
		 *(lambdaeps+0)=*(lambda+0)+eps;
		 *(lambdaeps+1)=0.0;
		 feps=ilyushin(shell, ii, lambdaeps, NULL);


		 /*UPDATE*/
		 if((*(feps+0)-*(f+0))==0.0)
		 {
		   FLAG++;
		   break;
		 }
		 else
		 {
		  *(dlambda+0) = -*(f+0)*eps / (*(feps+0)-*(f+0));
		 }

		 free(feps);

		 *(lambda+0) += *(dlambda+0);

		 /*NEW JUDGE*/
		 free(f);
		 f=ilyushin(shell, ii, lambda, NULL);


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
	   while(residual>tolerance && iteration<maxiteration)
	   {
		 /*SENSITIVITY*/
		 *(lambdaeps+0)=0.0;
		 *(lambdaeps+1)=*(lambda+1)+eps;
		 feps=ilyushin(shell, ii, lambdaeps, NULL);

		 /*UPDATE*/
		 if((*(feps+1)-*(f+1))==0.0)
		 {
		   FLAG++;
		   break;
		 }
		 else
		 {
		  *(dlambda+1) = -*(f+1)*eps / (*(feps+1)-*(f+1));
		 }
		 free(feps);

		 *(lambda+1) += *(dlambda+1);

		 /*NEW JUDGE*/
		 free(f);
		 f=ilyushin(shell, ii, lambda, NULL);

		 residual = abs(*(f+1));
		 iteration++;
	   }
	   *(lambda1 + 0) = *(lambda + 0);
	   *(lambda1 + 1) = *(lambda + 1);

	   *(f1 + 0) = *(f + 0);
	   *(f1 + 1) = *(f + 1);


	   if(FLAG==3)
	   {
		 errormessage("RETURN-MAPPING FAILED./n");
		 *(lambda + 0) = 0.0;
		 *(lambda + 1) = 0.0;
       }
	   else if( *(lambda01 + 0) > 0.0 && *(lambda01 + 1) > 0.0 && *(f01 + 0) <  tolerance && *(f01 + 1) <  tolerance)
	   {
		 *(lambda + 0) = *(lambda01 + 0);
		 *(lambda + 1) = *(lambda01 + 1);
	   }
	   else if( *(lambda0 + 0) > 0.0 && abs(*(f0 + 0)) < tolerance && *(f0 + 1) < 0.0)
	   {
		 *(lambda + 0) = *(lambda0 + 0);
		 *(lambda + 1) = *(lambda0 + 1);
	   }
	   else if( *(lambda1 + 1) > 0.0 && abs(*(f1 + 1)) < tolerance && *(f1 + 0) < 0.0)
	   {
		 *(lambda + 0) = *(lambda1 + 0);
		 *(lambda + 1) = *(lambda1 + 1);
	   }

   }
   free(f);
   f=ilyushin(shell, ii, lambda, 1);

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
   free(lambda);

   return;
}



/*ELASTO-PLASTICITY*/
double mises(struct oshell* shell, int ii, int jj, double lambda, int UPDATEFLAG)
{
  int i,j;
  struct gausspoint2* gp2;
  double* qstress,* ostress;
  double f;
  double E,t,poi;
  double yinit,y;
  double alpha;

  double v[3],p[3];/*[Q]^T[C][Q],[Q]^T[A][Q]*/
  double O[3],Oinv[3];
  double a[3];/*df/dσ*/


  char str[800];

  gp2 = &(shell->gp[ii].gp2[jj]);

  qstress=(double*)malloc(3 * sizeof(double));
  ostress=(double*)malloc(3 * sizeof(double));

  /*          | 1 -1  0 |  */
  /*[Q] = 1/√2| 1  1  0 |  */
  /*          | 0  0 √2 |  */

  /*[Q]diag[a,b,c][Q]^T  */
  /*    | a+b a-b   0 |  */
  /*=1/2| a-b a+b   0 |  */
  /*    |   0   0  2c |  */

  /*f={σ}^T[A]{σ}*/
  /*df/dσ=2[A]{σ}*/

  /*{σ}={σtry}-2dλ[C][A]{σ}                      */
  /*([I]+2dλ[C][A]){σ}=[W]^-1{σ}={σtry}          */
  /*{σ}=[W]{σtry}                                */

  /*[Q]^T{σ}=[Q]^T{σtry}-[Q]^T2dλ[C][A][Q][Q]^T{σ}   */
  /*[Q](I+[Q]^T2dλ[C][A][Q])[Q]^T{σ}=[W]^-1{σ}={σtry}*/

  /*[W]^-1=[Q][O][Q^T]                           */
  /*[W]=[Q][O]^-1[Q^T])                          */
  /*[O]=I+[Q]^T2dλ[C][A][Q]=I+2dλ[Q]^T[C][Q][Q]^T[A][Q]=diag{1 1 1}+2dλEdiag{v}diag{p}*/


  /*{σtry},dλ:KNOWN*/
  /*{σ}:UNKNOWN,MATRIX[W] CALCULATION NEEDED*/
  /*{σ}=([I]+2dλ[C][A])^-1{σtry}=[W]{σtry}*/

  E=shell->sect->E;
  t=shell->sect->area;
  poi=shell->sect->poi;
  yinit=shell->sect->yieldinit;


  v[0]=E/(1-poi);
  v[1]=E/(1+poi);
  v[2]=E/(2.0*(1+poi));

  p[0]=0.5;
  p[1]=1.5;
  p[2]=3.0;

  /*[O]*/
  for (i = 0; i < 3; i++)
  {
	O[i] = 1.0 + 2.0*lambda*v[i]*p[i];
  }
  /*Oinv:[O]^-1*/
  for (i = 0; i < 3; i++)
  {
	Oinv[i] =  1/O[i];
  }

  /*[Q]^T{σtry-αtry}*/
  *(qstress+0)=( gp2->stress[0] + gp2->stress[1] )/sqrt(2.0);
  *(qstress+1)=(-gp2->stress[0] + gp2->stress[1] )/sqrt(2.0);
  *(qstress+2)=  gp2->stress[2];


  /*[O]^-1*(diag[Q^T,Q^T]){σtry-αtry}*/
  *(ostress+0) = Oinv[0]**(qstress+0);
  *(ostress+1) = Oinv[1]**(qstress+1);
  *(ostress+2) = Oinv[2]**(qstress+2);

  /*f={σtry}^T[W]^T[A][W]{σtry}*/
  f = (p[0]**(ostress+0)**(ostress+0) + p[1]**(ostress+1)**(ostress+1) + p[2]**(ostress+2)**(ostress+2));

  alpha = gp2->alpha + 2.0*lambda*sqrt(f) ;
  y = yieldstress(shell->sect,alpha,NULL,NULL);

  f -= pow(y,2.0);

  if(UPDATEFLAG == 1)/*UPDATE PARAMS*/
  {
	  /*(diag[Q,Q])[O]^-1*(diag[Q^T,Q^T]){σtry-αtry}*/
	  gp2->stress[0]=(*(ostress+0)-*(ostress+1))/sqrt(2.0);
	  gp2->stress[1]=(*(ostress+0)+*(ostress+1))/sqrt(2.0);
	  gp2->stress[2]= *(ostress+2);

	  //gp2->backstress[0]=(*(ostress+0)-*(ostress+1))/sqrt(2.0);
	  //gp2->backstress[1]=(*(ostress+0)+*(ostress+1))/sqrt(2.0);
	  //gp2->backstress[2]= *(ostress+2);

	  gp2->alpha = alpha;
	  gp2->f = f;
	  gp2->lambda = lambda;

	  /*df/dσ*/
	  a[0] = 2.0 * gp2->stress[0] - gp2->stress[1];
	  a[1] = 2.0 * gp2->stress[1] - gp2->stress[0];
	  a[2] = 6.0 * gp2->stress[2];

	  /*dλ*df/dσ*/
	  for (i = 0; i < 3; i++)
	  {
		  gp2->estrain[i] -= lambda*a[i];
		  gp2->pstrain[i] += lambda*a[i];
	  }
  }
  free(qstress);
  free(ostress);

  return f;
}


double** shellCconsistentmises(struct oshell shell, int ii)
{
  int jj,i,j;
  int nnod,ngp,nstress;

  struct gausspoint* gp;
  struct gausspoint2* gp2;

  double E,t,poi;
  double v[3],p[3];/*[Q]^T[C][Q],[Q]^T[A][Q]*/
  double O[3], Oinv[3];
  double a[3];

  double ** C;/*3x3*/
  double ** consistentC;/*6x6*/
  double ** W;
  double ** H, * g;
  double Halpha, galpha;
  double y, dy, ddy;
  double yinit;

  double* N;
  double M;
  double beta;
  double tolerance = 1.0e-6;/*SAME AS RETURNMAP*/

  char str[800];

  nnod = shell.nnod;
  ngp = shell.ngp;
  nstress = shell.nstress;

  E=shell.sect->E;
  t=shell.sect->area;
  poi=shell.sect->poi;
  yinit = shell.sect->yieldinit;

  v[0]=E/(1-poi);
  v[1]=E/(1+poi);
  v[2]=E/(2.0*(1+poi));

  p[0]=0.5;
  p[1]=1.5;
  p[2]=3.0;

	consistentC = (double**)malloc(nstress * sizeof(double*));
	for (i = 0; i < nstress; i++)
	{
	  *(consistentC + i) = (double*)malloc(nstress * sizeof(double));
	  for (j = 0; j < nstress; j++)
	  {
		*(*(consistentC + i) + j) = 0.0;
	  }
	}

	W = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
	  *(W + i) = (double*)malloc(3 * sizeof(double));
	}
	g  = (double*)malloc(3 * sizeof(double*));

	gp = &(shell.gp[ii]);
	for(jj = 0; jj < gp->ngp2; jj++)/*FOR EACH INTEGRATION POINT*/
	{
	  gp2 = &(gp->gp2[jj]);

	  C=planeC(shell);

	  /*[O]*/
	  for (i = 0; i < 3; i++)
	  {
		O[i] = 1.0 + 2.0*(gp2->lambda)*v[i]*p[i];
	  }
	  /*Oinv:[O]^-1*/
	  for (i = 0; i < 3; i++)
	  {
		Oinv[i] =  1/O[i];
	  }

	  *(*(W+0)+0) = 0.5*(Oinv[0]+Oinv[1]);
	  *(*(W+0)+1) = 0.5*(Oinv[0]-Oinv[1]);
	  *(*(W+0)+2) = 0.0;
	  *(*(W+1)+0) = *(*(W+0)+1);
	  *(*(W+1)+1) = *(*(W+0)+0);
	  *(*(W+1)+2) = 0.0;
	  *(*(W+2)+0) = 0.0;
	  *(*(W+2)+1) = 0.0;
	  *(*(W+2)+2) = Oinv[2];

	　H = matrixmatrix(W,C,3);

	  /*{df/dσ}={d(σAσ)/dσ}=2[A]{σ}={g}*/
	  a[0] = 2.0 * gp2->stress[0] - gp2->stress[1];
	  a[1] = 2.0 * gp2->stress[1] - gp2->stress[0];
	  a[2] = 6.0 * gp2->stress[2];

	  *(g+0) = a[0];
	  *(g+1) = a[1];
	  *(g+2) = a[2];


	  /*df/dq=d/dq(-σy^2)=-2*σy*σy'(q)/=beta*/
	  y = yieldstress(shell.sect, gp2->alpha, &dy, &ddy);
	  if(dy!=0)
	  {
		Halpha = 1.0/(dy-(gp2->lambda)*2.0*dy*dy);
		galpha = -2.0*y*dy;

		beta = Halpha*pow(galpha,2.0);
	  }
	  else
	  {
		beta = 0;
	  }

	  if(abs(gp2->f)<tolerance)
	  {
		N=matrixvector(H,g,3);
		M=dotproduct(g,N,3);
		M+=beta;

		for (i = 0; i < 3; i++)
		{
		  for (j = 0; j < 3; j++)
		  {
			*(*(C+i)+j) = *(*(H+i)+j) - *(N+i)**(N+j)/M;
		  }
		}
		free(N);
	  }
	  else
	  {
		for (i = 0; i < 3; i++)
		{
		  for (j = 0; j < 3; j++)
		  {
			*(*(C+i)+j) = *(*(H+i)+j);
		  }
		}
	  }
	  freematrix(H,3);

	  /*[A B]*/
	  /*[B D]*/

	  /*[A]=Σ[C]*w*t/2*/
	  /*[B]=Σ[C]*z*w*t/2*t/2*/
	  /*[D]=Σ[C]*z*z*w*t/2*t/2*t/2*/
	  for (i = 0; i < 3; i++)
	  {
		for (j = 0; j < 3; j++)
		{
		  *(*(consistentC+  i)  +j) += *(*(C+i)+j)*(gp2->w)*t/2.0;
		  *(*(consistentC+3+i)  +j) += *(*(C+i)+j)*(gp2->z)*(gp2->w)*pow(t/2.0,2);
		  *(*(consistentC  +i)+3+j) += *(*(C+i)+j)*(gp2->z)*(gp2->w)*pow(t/2.0,2);
		  *(*(consistentC+3+i)+3+j) += *(*(C+i)+j)*(gp2->z)*(gp2->z)*(gp2->w)*pow(t/2.0,3);
		}
	  }
	  freematrix(C,3);
	}
	freematrix(W,3);
	free(g);

	return consistentC;
}




void returnmapmises(struct oshell* shell, int ii, int jj)
{
   int i,j;
   double lambda, lambdaeps, dlambda;
   double f, feps;
   double residual;
   double tolerance = 1.0e-6;/*SAME AS CONSISTENT_C*/
   double eps = 1.0e-2;
   char str[800];
   int iteration;
   int maxiteration=20;
   int FLAG=0;

   lambda = 0.0;
   f=mises(shell, ii, jj, lambda, NULL);/*MISES YIELD FUNCTION.{σ-α}^T[A]{σ-α}*/

   if(f>0)/*RETURN-MAPPING PROCEDURE USING NEWTON-RAPTHON METHOD FOR MISES YIELD CONDITION.*/
   {
	   /*YIELD SURFACE IS ACTIVE.*//*f=0 && lambda>0*/
	   residual=1.0;
	   iteration=1;
	   while(residual>tolerance && iteration<maxiteration)
	   {
		 /*SENSITIVITY*/
		 lambdaeps=lambda+eps;
		 feps=mises(shell, ii, jj, lambdaeps, NULL);

		 /*UPDATE*/
		 if(feps-f==0.0)
		 {
		   FLAG++;
		   break;
		 }
		 else
		 {
		   dlambda = -f*eps/(feps-f);
		 }
		 lambda += dlambda;

		 /*NEW JUDGE*/
		 f=mises(shell, ii, jj, lambda, NULL);

		 residual = abs(f);
		 iteration++;
	   }
	   if(FLAG==1)
	   {
		 errormessage("RETURN-MAPPING FAILED.");
		 lambda=0.0;
       }

   }
   f=mises(shell, ii, jj, lambda, 1);/*OUTPUT TO OSHELL.*/
   return;
}





void assemshellestrain(struct oshell* shell, double*** B, double* edisp)
{
  int ii,jj,i;
  int nnod,ngp,ngp2,nstress;
  double t;
  double* gpstrain;

  nnod = shell->nnod;
  ngp = shell->ngp;
  nstress = shell->nstress;
  t=shell->sect->area;

  for(ii = 0; ii < ngp; ii++)/*FOR EACH INTEGRATION POINT*/
  {
	ngp2 = (shell->gp[ii]).ngp2;
	gpstrain = matrixvectorIII(*(B+ii),edisp,nstress,6*nnod);
	for(i = 0; i < 6; i++)
	{
	  (shell->gp[ii]).estrain[i] = *(gpstrain+i);
	  /*STRAIN {εx εy εxy γx γy γxy} FOR INTEGRATION POINT ii*/
	}
	free(gpstrain);

	if(ngp2!=0)
	{
		for(jj = 0; jj < ngp2; jj++)
		{
			for(i = 0; i < 3; i++)
			{
			  (shell->gp[ii]).gp2[jj].estrain[i] = (shell->gp[ii]).estrain[i]
												 + (shell->gp[ii]).estrain[3+i] * (shell->gp[ii]).gp2[jj].z * t/2.0;
			  /*STRAIN {εx εy εxy} FOR INTEGRATION POINT ii*//*-h/2<z<h/2*/
			}
		}
	}
  }
  return;
}


void assemshellestress(struct oshell* shell)
{
  int ii,jj,i;
  int nnod,ngp,nstress;
  int ngp2;
  double t;
  struct gausspoint* gp;
  struct gausspoint2* gp2;
  double* gpstrain;
  double* gpstress;
  double** C;

  nnod = shell->nnod;
  ngp = shell->ngp;
  nstress = shell->nstress;
  t=shell->sect->area;


  for(ii = 0; ii < ngp; ii++)/*FOR EACH INTEGRATION POINT*/
  {
	gp = &(shell->gp[ii]);
	ngp2 = gp->ngp2;

	if(ngp2==0)/*ILYUSHIN*/
	{
		C=shellC(*shell);

		gpstrain = (double*)malloc(nstress * sizeof(double));
		for(i = 0; i < nstress; i++)
		{
		  /*estrain : NEW TOTAL STRAIN DERIVED FROM [B]{d}*/
		  /*pstrain : OLD PLASTIC STRAIN*/
		  gp->estrain[i] -= gp->pstrain[i];/*UPDATE estrain : NEW ELASTIC STRAIN*/

		  *(gpstrain+i) = gp->estrain[i];
		}
		gpstress = matrixvector(C,gpstrain,nstress);/*ELASTIC PREDICTOR*/
		for(i = 0; i < nstress; i++)
		{
		  gp->stress[i] = *(gpstress + i);/*TRIAL STRESS*/
		}
		free(gpstress);
		free(gpstrain);

		/*RETURN-MAPPING & UPDATE*/
		if(shell->sect->yieldinit != 0.0)returnmapilyushin(shell, ii);

		/*
		(shell->gp[ii]).Ee = 0.0;
		for(i = 0; i < 6; i++)
		{
		  (shell->gp[ii]).Ee += 0.5 * (shell->area) * (shell->gp[ii]).w
							   * (shell->gp[ii]). stress[i]
							   * (shell->gp[ii]).estrain[i];
		  (shell->gp[ii]).Ep += (shell->area) * (shell->gp[ii]).w
							   * 0.5 *((shell->gp[ii]). stress[i] + (mshell->gp[ii]). stress[i])
							   *((shell->gp[ii]).pstrain[i] - (mshell->gp[ii]).pstrain[i]);
		}
		*/

		freematrix(C,6);
	}
    else/*MISES*/
	{
		C=planeC(*shell);

		for(i = 0; i < nstress; i++)
		{
		  gp->stress[i]　= 0.0;
		}

		for(jj = 0; jj < ngp2; jj++)
		{
			gp2 = &(gp->gp2[jj]);

			gpstrain = (double*)malloc(3 * sizeof(double));

			for(i = 0; i < 3; i++)
			{
			  /*estrain : NEW TOTAL STRAIN DERIVED FROM [B]{d}*/
			  /*pstrain : OLD PLASTIC STRAIN*/
			  gp2->estrain[i] -= gp2->pstrain[i];/*UPDATE estrain : NEW ELASTIC STRAIN*/

			  *(gpstrain+i) = gp2->estrain[i];
			}
			gpstress = matrixvector(C,gpstrain,3);/*ELASTIC PREDICTOR*/
			for(i = 0; i < 3; i++)
			{
			  gp2->stress[i] = *(gpstress + i);/*TRIAL STRESS*/
			}
			free(gpstress);
			free(gpstrain);

			/*RETURN-MAPPING & UPDATE*/
			if(shell->sect->yieldinit != 0.0)returnmapmises(shell, ii, jj);

			/*
			(shell->gp[ii]).gp2[jj].Ee = 0.0;
			for(i = 0; i < 3; i++)
			{
			  (shell->gp[ii]).gp2[jj].Ee += 0.5 * (shell->area) * (shell->gp[ii]).gp2[jj].w
										   * (shell->gp[ii]).gp2[jj].stress[i]
										   * (shell->gp[ii]).gp2[jj].estrain[i];
			  (shell->gp[ii]).gp2[jj].Ep += (shell->area) * (shell->gp[ii]).gp2[jj].w
										   * 0.5 *((shell->gp[ii]). stress[i] + (mshell->gp[ii]). stress[i])
										   *((shell->gp[ii]).pstrain[i] - (mshell->gp[ii]).pstrain[i]);
			}
			*/
			for(i = 0; i < 3; i++)
			{
				gp->stress[i]   += gp2->stress[i] * gp2->w * t/2.0;
				gp->stress[3+i] += gp2->stress[i] * gp2->z * gp2->w * pow(t/2.0,2);
			}
		}




		freematrix(C,3);
	}
  }
  return;
}





double* assemshelleinternal(struct oshell* shell, double*** B)
{
  int ii,jj,i,j;
  int nnod,ngp,ngp2,nstress;
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
	  *(einternal+i) += *(gpinternal+i)*(shell->gp[ii]).w*(shell->area);
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








