/*FUNCTIONS FOR SHELL*/
double*** elasticCshell(struct oshell shell);
double** shelldrccos(struct oshell shell);
double shellarea(struct oshell shell);

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

/*FOR PLASTIC (ILYUSHIN'S STRESS RESULTANT)*/
double* ilyushin(struct oshell shell, double* lambda, double* estress_try, double* estress, double** W);/*YIELD FUNCTION*/
double* returnmapilyushin(struct oshell shell, double** C, double* estress);/*RETURN-MAPPING & UPDATE STIFFNESS MATRIX*/
void consistentCilyushin(struct oshell shell, double** C, double** W, double* estress, double* lambda);/*UPDATE STIFFNESS MATRIX*/

double** assemshellestrain(struct oshell shell, double*** B, double* edisp);
double** assemshellestress(struct oshell shell, double*** C, double** estrain);
double* assemshelleinternal(struct oshell shell, double*** B, double** estress);


double*** elasticCshell(struct oshell shell)
{
	int i,j,ii;
	double*** C;
	double E = shell.sect->E;
	double t = shell.sect->area;
	double poi = shell.sect->poi;

	C = (double***)malloc(7 * sizeof(double**));
	for (ii = 0; ii < 7; ii++)
	{
		*(C+ii) = (double**)malloc(6 * sizeof(double*));
		for (i = 0; i < 6; i++)
		{
			*(*(C+ii) + i) = (double*)malloc(6 * sizeof(double));
			for (j = 0; j < 6; j++)
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

double*** assemshellshape(struct oshell shell, double** drccos)
{
  int i,j,k,ii;
  double *a,*b,*c;
  double *Lx,*Ly;
  double **exy;
  double area;
  double** L;
  double *aa,*bb,*cc,*dd,*ee;
  double len;
  double** Bp,**Bb;
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
  L=(double **)malloc(7*sizeof(double *));/*TRIANGULAR COORDINATION OF INTEGRATION POINTS*/
  for(i=0;i<7;i++)
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
  *(*(L+5)+0)=0.0;*(*(L+5)+1)=1.0;*(*(L+5)+2)=0.0;
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
  area = 0.5 * ( *(*(exy+1)+0)**(*(exy+2)+1) - *(*(exy+2)+0)**(*(exy+1)+1) );

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

  B=(double ***)malloc(7*sizeof(double **));
  for(ii=0;ii<7;ii++)
  {
	*(B+ii)=(double **)malloc(6*sizeof(double*));
	for(i=0;i<6;i++)
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
  for(ii=0;ii<7;ii++)
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
  }


  free(a);
  free(b);
  free(c);
  free(Lx);
  free(Ly);
  freematrix(exy,3);
  freematrix(L,7);
  free(aa);
  free(bb);
  free(cc);
  free(dd);
  free(ee);
  freematrix(Bp,3);
  freematrix(Bb,3);

  return B;
}


double **assemshellemtx(struct oshell shell)
/*ASSEMBLAGE ELASTIC MATRIX.*/
{
  int i,j;
  int ii;/*LOOP FOR INTEGRATION POINTS*/
  double** Ke;
  double* w;
  double area;
  double** drccos;
  double*** C, ***B;
  double** Dp,** Bp,** Bpt,** DpBp,** Kp;
  double** Db,** Bb,** Bbt,** DbBb,** Kb;

  double prate = shell.prate;
  double brate = shell.brate;

  Ke=(double **)malloc(18*sizeof(double *));
  for(i=0;i<18;i++)
  {
	*(Ke+i)=(double *)malloc(18*sizeof(double));
	for(j=0;j<18;j++)
	{
	  *(*(Ke+i)+j)=0.0;                                              /*INITIAL.*/
	}
  }

  w=(double *)malloc(7*sizeof(double));

  Dp=(double **)malloc(3*sizeof(double *));
  Db=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Dp+i)=(double *)malloc(3*sizeof(double));
	*(Db+i)=(double *)malloc(3*sizeof(double));
  }

  Bp=(double **)malloc(3*sizeof(double *));
  Bb=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bp+i)=(double *)malloc(6*sizeof(double));
	*(Bb+i)=(double *)malloc(9*sizeof(double));
  }


  *(w+0)=27.0/60.0;
  *(w+1)=8.0/60.0;
  *(w+2)=*(w+1);
  *(w+3)=*(w+1);
  *(w+4)=3.0/60.0;
  *(w+5)=*(w+4);
  *(w+6)=*(w+4);

  area = shellarea(shell);

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

  /*IN=PLANE*/
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bp+i)+3*j+0) = *(*(*(B+0)+i)+6*j+0);
		*(*(Bp+i)+3*j+1) = *(*(*(B+0)+i)+6*j+1);
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

  /*BENDING*/
  for(ii=0;ii<7;ii++)
  {
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
		*(*(Ke+6*i+2)+6*j+2) += *(*(Kb+3*i+0)+3*j+0)**(w+ii)*area*brate;
		*(*(Ke+6*i+2)+6*j+3) += *(*(Kb+3*i+0)+3*j+1)**(w+ii)*area*brate;
		*(*(Ke+6*i+2)+6*j+4) += *(*(Kb+3*i+0)+3*j+2)**(w+ii)*area*brate;
		*(*(Ke+6*i+3)+6*j+2) += *(*(Kb+3*i+1)+3*j+0)**(w+ii)*area*brate;
		*(*(Ke+6*i+3)+6*j+3) += *(*(Kb+3*i+1)+3*j+1)**(w+ii)*area*brate;
		*(*(Ke+6*i+3)+6*j+4) += *(*(Kb+3*i+1)+3*j+2)**(w+ii)*area*brate;
		*(*(Ke+6*i+4)+6*j+2) += *(*(Kb+3*i+2)+3*j+0)**(w+ii)*area*brate;
		*(*(Ke+6*i+4)+6*j+3) += *(*(Kb+3*i+2)+3*j+1)**(w+ii)*area*brate;
		*(*(Ke+6*i+4)+6*j+4) += *(*(Kb+3*i+2)+3*j+2)**(w+ii)*area*brate;
	  }
	}

	freematrix(Bbt,9);
	freematrix(DbBb,3);
	freematrix(Kb,9);
  }

  /*DRILL*/
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

  free(w);
  freematrix(drccos,3);
  for(ii=0;ii<7;ii++)freematrix(*(C+ii), 6);
  free(C);
  for(ii=0;ii<7;ii++)freematrix(*(B+ii), 6);
  free(B);

  freematrix(Dp,3);
  freematrix(Db,3);

  freematrix(Bp,3);
  freematrix(Bb,3);

  return Ke;
}


double **assemshellpmtx(struct oshell shell,double*** C,double *** B)
/*ASSEMBLAGE ELASTO-PLASTIC MATRIX.*/
{
  int i,j;
  int ii;/*LOOP FOR INTEGRATION POINTS*/
  double** K;
  double* w;
  double area;
  double** Dp,** Bp,** Bpt,** DpBp,** Kp;
  double** Db,** Bb,** Bbt,** DbBb,** Kb;

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

  w=(double *)malloc(7*sizeof(double));

  Dp=(double **)malloc(3*sizeof(double *));
  Db=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Dp+i)=(double *)malloc(3*sizeof(double));
	*(Db+i)=(double *)malloc(3*sizeof(double));
  }

  Bp=(double **)malloc(3*sizeof(double *));
  Bb=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bp+i)=(double *)malloc(6*sizeof(double));
	*(Bb+i)=(double *)malloc(9*sizeof(double));
  }


  *(w+0)=27.0/60.0;
  *(w+1)=8.0/60.0;
  *(w+2)=*(w+1);
  *(w+3)=*(w+1);
  *(w+4)=3.0/60.0;
  *(w+5)=*(w+4);
  *(w+6)=*(w+4);

  area = shellarea(shell);


  for(ii=0;ii<7;ii++)
  {
	/*STRESS-STRAIN MATRIX*/
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Dp+i)+j) = *(*(*(C+ii)+i)+j);
		*(*(Db+i)+j) = *(*(*(C+ii)+3+i)+3+j);
	  }
	}

	/*IN=PLANE*/
	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(Bp+i)+3*j+0) = *(*(*(B+ii)+i)+6*j+0);
		*(*(Bp+i)+3*j+1) = *(*(*(B+ii)+i)+6*j+1);
	  }
	}

	Bpt=matrixtransposeIII(Bp,3,6);
	DpBp=matrixmatrixIII(Dp,Bp,3,3,6);
	Kp=matrixmatrixIII(Bpt,DpBp,6,3,6);


	for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	  {
		*(*(K+6*i+0)+6*j+0)+=*(*(Kp+2*i+0)+2*j+0)**(w+ii)*area*prate;
		*(*(K+6*i+0)+6*j+1)+=*(*(Kp+2*i+0)+2*j+1)**(w+ii)*area*prate;
		*(*(K+6*i+1)+6*j+0)+=*(*(Kp+2*i+1)+2*j+0)**(w+ii)*area*prate;
		*(*(K+6*i+1)+6*j+1)+=*(*(Kp+2*i+1)+2*j+1)**(w+ii)*area*prate;
	  }
	}

	freematrix(Bpt,6);
	freematrix(DpBp,3);
	freematrix(Kp,6);

	/*BENDING*/
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
		*(*(K+6*i+2)+6*j+2) += *(*(Kb+3*i+0)+3*j+0)**(w+ii)*area*brate;
		*(*(K+6*i+2)+6*j+3) += *(*(Kb+3*i+0)+3*j+1)**(w+ii)*area*brate;
		*(*(K+6*i+2)+6*j+4) += *(*(Kb+3*i+0)+3*j+2)**(w+ii)*area*brate;
		*(*(K+6*i+3)+6*j+2) += *(*(Kb+3*i+1)+3*j+0)**(w+ii)*area*brate;
		*(*(K+6*i+3)+6*j+3) += *(*(Kb+3*i+1)+3*j+1)**(w+ii)*area*brate;
		*(*(K+6*i+3)+6*j+4) += *(*(Kb+3*i+1)+3*j+2)**(w+ii)*area*brate;
		*(*(K+6*i+4)+6*j+2) += *(*(Kb+3*i+2)+3*j+0)**(w+ii)*area*brate;
		*(*(K+6*i+4)+6*j+3) += *(*(Kb+3*i+2)+3*j+1)**(w+ii)*area*brate;
		*(*(K+6*i+4)+6*j+4) += *(*(Kb+3*i+2)+3*j+2)**(w+ii)*area*brate;
	  }
	}

	freematrix(Bbt,9);
	freematrix(DbBb,3);
	freematrix(Kb,9);
  }

  /*DRILL*/
  double E=shell.sect->E;
  double t=shell.sect->area;

  for(ii=0;ii<1;ii++)
  {
	*(*(K+ 5)+ 5)= 0.030*E*t*area;
	*(*(K+ 5)+11)=-0.015*E*t*area;
	*(*(K+ 5)+17)=*(*(K+5)+11);
	*(*(K+11)+ 5)=*(*(K+5)+11);
	*(*(K+11)+11)=*(*(K+5)+ 5);
	*(*(K+11)+17)=*(*(K+5)+11);
	*(*(K+17)+ 5)=*(*(K+5)+11);
	*(*(K+17)+11)=*(*(K+5)+11);
	*(*(K+17)+17)=*(*(K+5)+ 5);
  }

  free(w);
  freematrix(Dp,3);
  freematrix(Db,3);
  freematrix(Bp,3);
  freematrix(Bb,3);

  return K;
}



/*
void addkpemtx(double **Ke,struct oshell shell,double** C,double ***B)
{
  int i,j;
  double **Dp,**Bp,**Bpt,**DpBp,**Kp;
  double A;
  double prate = shell.prate;

  Dp=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Dp+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	  *(*(Dp+i)+j)=*(*(C+i)+j);
	}
  }
  Bp=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bp+i)=(double *)malloc(6*sizeof(double));
	for(j=0;j<3;j++)
	{
	  *(*(Bp+i)+3*j+0) = *(*(*(B+ii)+i)+6*j+0);
	  *(*(Bp+i)+3*j+1) = *(*(*(B+ii)+i)+6*j+1);
	}
  }

  Bpt=matrixtransposeIII(Bp,3,6);
  DpBp=matrixmatrixIII(Dp,Bp,3,3,6);
  Kp=matrixmatrixIII(Bpt,DpBp,6,3,6);


  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Ke+6*i+0)+6*j+0)=*(*(Kp+2*i+0)+2*j+0)*A*prate;
	  *(*(Ke+6*i+0)+6*j+1)=*(*(Kp+2*i+0)+2*j+1)*A*prate;
	  *(*(Ke+6*i+1)+6*j+0)=*(*(Kp+2*i+1)+2*j+0)*A*prate;
	  *(*(Ke+6*i+1)+6*j+1)=*(*(Kp+2*i+1)+2*j+1)*A*prate;
	}
  }

  freematrix(Dp,3);
  freematrix(Bp,3);
  freematrix(Bpt,6);
  freematrix(DpBp,3);
  freematrix(Kp,6);
  return;
}

void addkbemtx(double **Ke,struct oshell shell,double** C,double ***B)
{
  int i,j,ii;
  double* w;
  double **Db,**Bb,**Bbt,**DbBb,**Kb;
  double A;
  double brate = shell.brate;

  w=(double *)malloc(7*sizeof(double));

  *(w+0)=27.0/60.0;
  *(w+1)=8.0/60.0;
  *(w+2)=8.0/60.0;
  *(w+3)=8.0/60.0;
  *(w+4)=3.0/60.0;
  *(w+5)=3.0/60.0;
  *(w+6)=3.0/60.0;

  Db=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Db+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	  *(*(Db+i)+j) = *(*(C+3+i)+3+j)
	}
  }
  Bb=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Bb+i)=(double *)malloc(9*sizeof(double));
  }


  for(ii=0;ii<7;ii++)
  {
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
		*(*(Ke+6*i+2)+6*j+2) += *(*(Kb+3*i+0)+3*j+0)**(w+ii)*A*brate;
		*(*(Ke+6*i+2)+6*j+3) += *(*(Kb+3*i+0)+3*j+1)**(w+ii)*A*brate;
		*(*(Ke+6*i+2)+6*j+4) += *(*(Kb+3*i+0)+3*j+2)**(w+ii)*A*brate;
		*(*(Ke+6*i+3)+6*j+2) += *(*(Kb+3*i+1)+3*j+0)**(w+ii)*A*brate;
		*(*(Ke+6*i+3)+6*j+3) += *(*(Kb+3*i+1)+3*j+1)**(w+ii)*A*brate;
		*(*(Ke+6*i+3)+6*j+4) += *(*(Kb+3*i+1)+3*j+2)**(w+ii)*A*brate;
		*(*(Ke+6*i+4)+6*j+2) += *(*(Kb+3*i+2)+3*j+0)**(w+ii)*A*brate;
		*(*(Ke+6*i+4)+6*j+3) += *(*(Kb+3*i+2)+3*j+1)**(w+ii)*A*brate;
		*(*(Ke+6*i+4)+6*j+4) += *(*(Kb+3*i+2)+3*j+2)**(w+ii)*A*brate;
	  }
	}


	freematrix(Bbt,9);
	freematrix(DbBb,3);
	freematrix(Kb,9);
  }
  freematrix(Bb,3);
  freematrix(Db,3);
  free(w);

  return;
}

void addktemtx(double **Ke,struct oshell shell)
{
	double E,t;

	E=shell.sect->E;
	t=shell.sect->area;

	*(*(Ke+5)+5)=0.03*E*t*A;
	*(*(Ke+5)+11)=-0.015*E*t*A;
	*(*(Ke+5)+17)=*(*(Ke+5)+11);
	*(*(Ke+11)+5)=*(*(Ke+5)+11);
	*(*(Ke+11)+11)=*(*(Ke+5)+5);
	*(*(Ke+11)+17)=*(*(Ke+5)+11);
	*(*(Ke+17)+5)=*(*(Ke+5)+11);
	*(*(Ke+17)+11)=*(*(Ke+5)+11);
	*(*(Ke+17)+17)=*(*(Ke+5)+5);
	return;
}
*/





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



double* ilyushin(struct oshell shell, double* lambda, double* estress_try, double* estress, double** W)
{
  int i,j;
  double* qstress,* ostress;
  double* f;
  double fn,fm;
  double fn2,fm2,fnfm;
  double E,t,A,I,v;
  double lambda1,lambda2;

  double c[3],p[3],e[3],det[3];
  double O[2][2][3], Oinv[2][2][3];

  qstress=(double*)malloc(6 * sizeof(double));
  ostress=(double*)malloc(6 * sizeof(double));
  f=(double*)malloc(2 * sizeof(double));

  if(lambda==NULL)
  {
	lambda1=0.0;
	lambda2=0.0;
  }
  else
  {
	lambda1=*(lambda+0);
	lambda2=*(lambda+1);
  }


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

  E=shell.sect->E;
  t=shell.sect->area;
  v=shell.sect->poi;

  fn=shell.sect->fmax[0];/*fy*t*/
  fm=shell.sect->fmax[3];/*fy*pow(t,3)/4.0;*/

  fn2=pow(fn,2);
  fm2=pow(fm,2);
  fnfm=fn*fm;

  A=E*t;
  I=E*pow(t,3)/12.0;

  c[0]=1/(1-v);
  c[1]=1/(1+v);
  c[2]=1/(2.0*(1+v));

  p[0]=0.5;
  p[1]=1.5;
  p[2]=3.0;

  e[0] = (lambda1 + lambda2) / fn2;
  e[1] = (lambda1 + lambda2) / fm2;
  e[2] = (lambda1 - lambda2) /(2.0*sqrt(3.0)*fnfm);

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

  /*(diag[Q^T,Q^T]){ƒÐtry-ƒ¿try}*/
  *(qstress+0)=( *(estress_try+0)+*(estress_try+1))/sqrt(2.0);
  *(qstress+1)=(-*(estress_try+0)+*(estress_try+1))/sqrt(2.0);
  *(qstress+2)=  *(estress_try+2);
  *(qstress+3)=( *(estress_try+3)+*(estress_try+4))/sqrt(2.0);
  *(qstress+4)=(-*(estress_try+3)+*(estress_try+4))/sqrt(2.0);
  *(qstress+5)=  *(estress_try+5);

  /*[O]^-1*(diag[Q^T,Q^T]){ƒÐtry-ƒ¿try}*/
  *(ostress+0) = Oinv[0][0][0]**(qstress+0) + Oinv[0][1][0]**(qstress+3);
  *(ostress+1) = Oinv[0][0][1]**(qstress+1) + Oinv[0][1][1]**(qstress+4);
  *(ostress+2) = Oinv[0][0][2]**(qstress+2) + Oinv[0][1][2]**(qstress+5);
  *(ostress+3) = Oinv[1][0][0]**(qstress+0) + Oinv[1][1][0]**(qstress+3);
  *(ostress+4) = Oinv[1][0][1]**(qstress+1) + Oinv[1][1][0]**(qstress+4);
  *(ostress+5) = Oinv[1][0][2]**(qstress+2) + Oinv[1][1][0]**(qstress+5);

  if(estress!=NULL)
  {
	/*(diag[Q,Q])[O]^-1*(diag[Q^T,Q^T]){ƒÐtry-ƒ¿try}*/
	*(estress+0)=(*(ostress+0)-*(ostress+1))/sqrt(2.0);
	*(estress+1)=(*(ostress+0)+*(ostress+1))/sqrt(2.0);
	*(estress+2)= *(ostress+2);
	*(estress+3)=(*(ostress+3)-*(ostress+4))/sqrt(2.0);
	*(estress+4)=(*(ostress+3)+*(ostress+4))/sqrt(2.0);
	*(estress+5)= *(ostress+5);
  }

  if(W!=NULL)
  {
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
  }

  q = shell.q;
  yinit = shell.yinit;

  y = isohardening(q);



  /*f={ƒÐtry-ƒ¿try}^T[W]^T[A][W]{ƒÐtry-ƒ¿try}*/
  *(f+0) = p[0] * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 + *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + p[1] * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 + *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + p[2] * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 + *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm))
		 - pow(y/yinit,2);
  *(f+1) = p[0] * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 - *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + p[1] * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 - *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + p[2] * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 - *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm))
		 - pow(y/yinit,2);

  return f;
}


double* returnmapilyushin(struct oshell shell, double** C, double* estress)
{
   int i,j;
   double* lambda,* lambdaeps,* dlambda;
   double* f,* feps;
   double* estress_try;
   double** dfdl;
   double** W;
   double det;
   double residual;
   double tolerance = 1.0e-8;
   double eps = 1.0e-8;
   int FLAG;

   /*FLAG = 0 ; BOTH SURFACES ARE ACTIVE*/
   /*FLAG = 1 ; SURFACE 1 IS ACTIVE     */
   /*FLAG = 2 ; SURFACE 2 IS ACTIVE     */


   lambda = (double*)malloc(2 * sizeof(double));

   /*YIELD JUDGEMENT.*/
   f=ilyushin(shell, NULL, estress, NULL, NULL);/*ILYUSHIN'S YIELD FUNCTION.{ƒÐ-ƒ¿}^T[A]{ƒÐ-ƒ¿}*/

   if( *(f+0)>0 || *(f+1)>0 )/*IF YIELD SURFACE ACTIVE.*/
   {
	   dfdl= (double**)malloc(2 * sizeof(double*));
	   for(i=0;i<2;i++)
	   {
		 *(dfdl+i)= (double*)malloc(2 * sizeof(double));
	   }
	   dlambda = (double*)malloc(2 * sizeof(double));
	   lambdaeps = (double*)malloc(2 * sizeof(double));
	   estress_try = (double*)malloc(6 * sizeof(double));

	   /*RETURN-MAPPING PROCEDURE USING NEWTON-RAPTHON METHOD FOR ILYUSHIN'S YIELD CONDITION.*/
		for(i=0;i<2;i++)
		{
		  *(lambda+i)=0.0;
		}
		for(i=0;i<6;i++)
		{
		  *(estress_try+i)=*(estress+i);
		}

		FLAG = 0;

		while(residual>tolerance)/*BOTH YIELD SURFACES ARE ACTIVE.*/
		{
		  residual = 0.0;
		  /*YIELD FUNCTION*/
		  f=ilyushin(shell, lambda, estress_try, NULL, NULL);
		  /*SENSITIVITY*/
		  for(i=0;i<2;i++)
		  {
			for(j=0;j<2;j++)
			{
			  if(i==j) *(lambdaeps+i)=*(lambda+i)+eps;
			  else     *(lambdaeps+i)=*(lambda+i);
			}
			feps=ilyushin(shell, lambdaeps, estress_try, NULL, NULL);
			for(j=0;j<2;j++)
			{
			  *(*(dfdl+j)+i)=(*(feps+j)-*(f+j))/eps;
			}
			free(feps);
		  }
		  free(f);

		  det=*(*(dfdl+0)+0)**(*(dfdl+1)+1)-*(*(dfdl+0)+1)**(*(dfdl+1)+0);

		  *(dlambda + 0) = -( *(*(dfdl+1)+1)**(f+0)-*(*(dfdl+0)+1)**(f+1))/det;
		  *(dlambda + 1) = -(-*(*(dfdl+1)+0)**(f+0)+*(*(dfdl+0)+0)**(f+1))/det;

		  if(*(dlambda + 0) < 0)
		  {
			  FLAG = 2;
			  break;
		  }
		  if(*(dlambda + 1) < 0)
		  {
			  FLAG = 1;
			  break;
		  }
		   *(lambda + 0) += *(dlambda + 0);
		   *(lambda + 1) += *(dlambda + 1);
		}
		if(FLAG == 1)
		{
			while(residual>tolerance)/*YIELD SURFACE 1 IS ACTIVE.*//*f1=0 && f2<0 && lambda1>0 && lambda2=0*/
			{
			  residual = 0.0;
			  /*YIELD FUNCTION*/
			  f=ilyushin(shell, lambda, estress_try, NULL, NULL);
			  /*SENSITIVITY*/
			  *(lambdaeps+0)=*(lambda+0)+eps;
			  *(lambdaeps+1)=0.0;

			  feps=ilyushin(shell, lambdaeps, estress_try, NULL, NULL);
			  *(feps+0)=(*(feps+0)-*(f+0))/eps;

			  *(lambda+0) += -*(f+0) / *(feps+0);

			  free(f);
			  free(feps);
			}
		}
		if(FLAG == 2)
		{
			while(residual>tolerance)/*YIELD SURFACE 2 IS ACTIVE.*//*f1=0 && f2<0 && lambda1>0 && lambda2=0*/
			{
			  residual = 0.0;
			  /*YIELD FUNCTION*/
			  f=ilyushin(shell, lambda, estress_try, NULL, NULL);

			  /*SENSITIVITY*/
			  *(lambdaeps+1)=*(lambda+1)+eps;
			  *(lambdaeps+0)=0.0;
			  feps=ilyushin(shell, lambdaeps, estress_try, NULL, NULL);
			  *(feps+1)=(*(feps+1)-*(f+1))/eps;

			  *(lambda+1) += -*(f+1) / *(feps+1);

			  free(f);
			  free(feps);
			}
		}

		/*UPDATE STRESS*/
		W = (double**)malloc(6 * sizeof(double*));
		for (i = 0; i < 6; i++)
		{
		  *(W + i) = (double*)malloc(6 * sizeof(double));
		}

		f = ilyushin(shell, lambda, estress_try, estress, W);

		consistentCilyushin(shell, mshell, C, W, estress, lambda);

		free(estress_try);
		free(lambdaeps);
		freematrix(dfdl,2);
   }
   else/*ELASTIC.*/
   {
		for(i=0;i<2;i++)
		{
			*(lambda+i)=0.0;
		}
   }
	return lambda;
}


void consistentCilyushin(struct oshell shell, struct memoryshell mshell, double** C, double** W, double* estress, double* lambda)
{
  int i,j;
  double fn,fm;
  double e[3],a[6];
  double** H;
  double* g1,* g2,* N1,* N2;
  double M11,M12,M21,M22;
  double det;

  fn=shell.sect->fmax[0];
  fn=shell.sect->fmax[3];

  H = matrixmatrix(W,C,6);
  for (i = 0; i < 6; i++)
  {
	for (j = 0; j < 6; j++)
	{
		*(*(C+i)+j) = *(*(H+i)+j);
	}
  }


  e[0] = fn*fn;
  e[1] = fm*fm;
  e[2] = 2.0*sqrt(3.0)*fn*fm;

  /*2[A]*/
  a[0] = 2.0**(estress+0) - *(estress+1);
  a[1] = 2.0**(estress+1) - *(estress+0);
  a[2] = 6.0**(estress+2);
  a[3] = 2.0**(estress+3) - *(estress+4);
  a[4] = 2.0**(estress+4) - *(estress+3);
  a[5] = 6.0**(estress+5);

  if(*(lambda+0)>0 && *(lambda+1)>0)
  {
	g1 = (double*)malloc(6 * sizeof(double));
	g2 = (double*)malloc(6 * sizeof(double));

	*(g1+0) = a[0]/e[0] + a[3]/e[2];
	*(g1+0) = a[1]/e[0] + a[4]/e[2];
	*(g1+0) = a[2]/e[0] + a[5]/e[2];
	*(g1+0) = a[0]/e[2] + a[3]/e[1];
	*(g1+0) = a[1]/e[2] + a[4]/e[1];
	*(g1+0) = a[2]/e[2] + a[5]/e[1];

	*(g2+0) =  a[0]/e[0] - a[3]/e[2];
	*(g2+0) =  a[1]/e[0] - a[4]/e[2];
	*(g2+0) =  a[2]/e[0] - a[5]/e[2];
	*(g2+0) = -a[0]/e[2] + a[3]/e[1];
	*(g2+0) = -a[1]/e[2] + a[4]/e[1];
	*(g2+0) = -a[2]/e[2] + a[5]/e[1];

	for (i = 0; i < 6; i++)
	{
	  *(g+i) = *(lambda+0)**(g1+i) + *(lambda+1)**(g2+i);
	}

	N1=matrixvector(H,g1,6);
	N2=matrixvector(H,g2,6);

	M11=dotproduct(g1,N1,6);
	M12=dotproduct(g1,N2,6);
	M21=M12;
	M22=dotproduct(g2,N2,6);

	det = M11*M22-M12*M21;

	for (i = 0; i < 6; i++)
	{
	  for (j = 0; j < 6; j++)
	  {
		*(*(C+i)+j) -= ( M22**(N1+i)**(N1+j) - M12**(N1+i)**(N2+j) - M21**(N2+i)**(N1+j) + M11**(N2+i)**(N2+j) )/det;
	  }
	}

	free(g1);
	free(g2);
	free(N1);
	free(N2);
  }
  else if(*(lambda+0)>0)
  {
	g1 = (double*)malloc(6 * sizeof(double));
	*(g1+0) = a[0]/e[0] + a[3]/e[2];
	*(g1+0) = a[1]/e[0] + a[4]/e[2];
	*(g1+0) = a[2]/e[0] + a[5]/e[2];
	*(g1+0) = a[0]/e[2] + a[3]/e[1];
	*(g1+0) = a[1]/e[2] + a[4]/e[1];
	*(g1+0) = a[2]/e[2] + a[5]/e[1];

	for (i = 0; i < 6; i++)
	{
	  *(g+i) = *(lambda+0)**(g1+i);
	}

	N1=matrixvector(H,g1,6);
	M11=dotproduct(g1,N1,6);

	for (i = 0; i < 6; i++)
	{
	  for (j = 0; j < 6; j++)
	  {
		*(*(C + i) + j) -= *(N1+i)**(N1+j)/M11;
	  }
	}

	free(g1);
	free(N1);
  }
  else if(*(lambda+1)>0)
  {
	g2 = (double*)malloc(6 * sizeof(double));
	*(g2+0) =  a[0]/e[0] - a[3]/e[2];
	*(g2+0) =  a[1]/e[0] - a[4]/e[2];
	*(g2+0) =  a[2]/e[0] - a[5]/e[2];
	*(g2+0) = -a[0]/e[2] + a[3]/e[1];
	*(g2+0) = -a[1]/e[2] + a[4]/e[1];
	*(g2+0) = -a[2]/e[2] + a[5]/e[1];

	for (i = 0; i < 6; i++)
	{
	  *(g+i) = *(lambda+1)**(g2+i);
	}

	N2=matrixvector(H,g2,6);
	M22=dotproduct(g2,N2,6);

	for (i = 0; i < 6; i++)
	{
	  for (j = 0; j < 6; j++)
	  {
		*(*(C + i) + j) -= *(N2+i)**(N2+j)/M22;
	  }
	}

	free(g2);
	free(N2);
  }

  return;
}




double** assemshellestrain(struct oshell shell, double*** B, double* edisp)
{
  int ii,i;
  int nnod;
  double* e;
  double** estrain;

  nnod = shell.nnod;
  estrain = (double**)malloc(7 * sizeof(double*));
  for(ii = 0; ii < 7; ii++)/*FOR EACH INTEGRATION POINT*/
  {
	 *(estrain+ii) = (double*)malloc(6 * sizeof(double));
  }

  for(ii = 0; ii < 7; ii++)/*FOR EACH INTEGRATION POINT*/
  {
	e = matrixvectorIII(*(B+ii),edisp,6,6*nnod); /*INCREMENTAL STRAIN d{ƒÃx ƒÃy ƒÃxy ƒÁx ƒÁy ƒÁxy} FOR NODE-ii*/
	for(i = 0; i < 6; i++)
	{
	  *(*(estrain+ii)+i) = *(e+i);
	}
	free(e);
  }
  return estrain;
}

double** assemshellestress(struct oshell shell, double*** C, double** estrain)
{
  int ii,i;
  int nnod;
  double* e;
  double** estress;
  double* lambda;

  nnod = shell.nnod;
  estress = (double**)malloc(7 * sizeof(double*));
  for(ii = 0; ii < 7; ii++)/*FOR EACH INTEGRATION POINT*/
  {
	 *(estress+ii) = (double*)malloc(6 * sizeof(double));
  }

  for(ii = 0; ii < 7; ii++)/*FOR EACH INTEGRATION POINT*/
  {
	for(i = 0; i < 6; i++)
	{
	  //*(*(estrain+ii)+i) -= shell->pstrain[ii][i];/*TRIAL STRESS*/
	}

	e = matrixvector(*(C+ii),*(estrain+ii),6);/*ELASTIC PREDICTOR*/


	/*UPDATE STRESS RESULTANT*/
	lambda = returnmapilyushin(shell, *(C+ii), e);

	for(i = 0; i < 6; i++)
	{
	  *(*(estress+ii)+i) = *(e+i);
	}

	free(e);
	free(lambda);
  }
  return estress;
}

double* assemshelleinternal(struct oshell shell, double*** B, double** estress)
{
  int ii,i;
  int nnod;
  double area;
  double* w;
  double* einternal,*e;
  double** Bt;

  nnod = shell.nnod;
  einternal = (double*)malloc(6 * nnod * sizeof(double));
  w = (double*)malloc(7 * sizeof(double));

  area = shellarea(shell);

  *(w+0)=27.0/60.0*area;
  *(w+1)=8.0/60.0*area;
  *(w+2)=*(w+1);
  *(w+3)=*(w+1);
  *(w+4)=3.0/60.0*area;
  *(w+5)=*(w+4);
  *(w+6)=*(w+4);

  for(ii = 0; ii < 7; ii++)
  {
	Bt = matrixtransposeIII(*(B+ii),6,6*nnod);
	e = matrixvectorIII(Bt,*(estress+ii),6*nnod,6);

	for(i = 0; i < 6; i++)
	{
	  *(einternal+i) += *(e+i)**(w+ii);
	}

	freematrix(Bt,6*nnod);
	free(e);
  }

  return einternal;
}

double isohardening(double q)/*equivalent plastic strain*/
{
  double y,yinit,c,n;
  yinit =  100;
  c = 100;
  n = 10;

  y = yinit + c * q;
  //y = yinit + c * pow(q,n);

  return y;
}

double gradisohardening(double q)
{
  double k,y,yinit,c,n;
  yinit =  100;
  c = 100;
  n = 10;

  k = c;
  //k = n * c * pow(q,(n-1));

  return k;
}
