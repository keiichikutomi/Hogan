/*FUNCTIONS FOR SHELL*/
double** elasticCshell(struct oshell shell);
double*** assemshellshapef(struct oshell shell, double** drccos);
double** assemshellemtx(struct oshell shell,double** C,double ***B);


double** elasticCshell(struct oshell shell)
{
	int i,j;
	double** C

	C = (double**)malloc(6 * sizeof(double*));
	for (i = 0; i < 6; i++)
	{
		*(C + i) = (double*)malloc(6 * sizeof(double));
		for (j = 0; j < 6; j++)
		{
			*(*(C + i) + j) = 0.0;
		}
	}

	*(*(C+0)+0)=E*t/(1-poi*poi);
	*(*(C+0)+1)=poi**(*(C+0)+0));
	*(*(C+1)+0)=*(*(C+0)+1);
	*(*(C+1)+1)=*(*(C+0)+0);
	*(*(C+2)+2)=0.5*E*t/(1+poi);

	*(*(C+3)+3)=E*pow(t,3)/(12.0*(1-poi*poi));
	*(*(C+3)+4)=poi**(*(C+0)+0));
	*(*(C+4)+3)=*(*(C+3)+4);
	*(*(C+4)+4)=*(*(C+3)+3);
	*(*(C+5)+5)=0.5*E*pow(t,3)/(12.0*(1+poi));

	return C;
}

double*** assemshellshapef(struct oshell shell, double** drccos)
{
  int i,j,k,ii;
  double *a,*b,*c;
  double *Lx,*Ly;
  double **exy;
  double** L;
  double *aa,*bb,*cc,*dd,*ee;
  double len;

  double** Bp,**Bb;


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

  for(i=0;i<3;i++)
  {
	  j=i+1;
	  k=i+2;
	  if(j>=3)j-=3;
	  if(k>=3)k-=3;

	  *(a+i)=*(*(exy+j)+0)**(*(exy+k)+1)-*(*(exy+k)+0)**(*(exy+j)+1);/*xjyk-xkyj*/
	  *(b+i)=*(*(exy+j)+1)-*(*(exy+k)+1);/*yj-yk*/
	  *(c+i)=*(*(exy+k)+0)-*(*(exy+j)+0);/*xk-xj*/

	  *(Lx+i)=0.5**(b+i)/A;/*dLi/dx=bi/(2*area)*/
	  *(Ly+i)=0.5**(c+i)/A;/*dLi/dy=ci/(2*area)*/
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
	/*Zienkiewiczらの非適合三角形要素T-9N*/
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
	/*Stricklin/Dhattらの離散Kirchhoff仮定三角形(DKT:Discrete Kirchhoff Triangular)要素T-9D*/
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

double **assemshellemtx(struct oshell shell,double** C,double ***B)
/*ASSEMBLAGE ELASTIC MATRIX.*/
{
  int i,j,ii;
  double **Ke

  Ke=(double **)malloc(18*sizeof(double *));
  for(i=0;i<18;i++)
  {
	*(Ke+i)=(double *)malloc(18*sizeof(double));
	for(j=0;j<18;j++)
	{
	  *(*(Ke+i)+j)=0.0;                                              /*INITIAL.*/
	}
  }

  //addkpemtx(Ke,shell,C,B);
  //addkbemtx(Ke,shell,C,B);
  //addktemtx(Ke,shell);


  double* w;
  double**Dp,**Bp,**Bpt,**DpBp,**Kp;
  double**Db,**Bb,**Bbt,**DbBb,**Kb;
  double E=shell.sect->E;
  double t=shell.sect->t;
  double A=shell.sect->area;
  double prate = shell.prate;
  double brate = shell.brate;

  w=(double *)malloc(7*sizeof(double));

  *(w+0)=27.0/60.0;
  *(w+1)=8.0/60.0;
  *(w+2)=8.0/60.0;
  *(w+3)=8.0/60.0;
  *(w+4)=3.0/60.0;
  *(w+5)=3.0/60.0;
  *(w+6)=3.0/60.0;

  Dp=(double **)malloc(3*sizeof(double *));
  Db=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Dp+i)=(double *)malloc(3*sizeof(double));
	*(Db+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	  *(*(Dp+i)+j) = *(*(C+i)+j);
	  *(*(Db+i)+j) = *(*(C+3+i)+3+j)
	}
  }

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

  /*IN=PLANE*/
  for(ii=0;ii<1;ii++)
  {
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
		*(*(Ke+6*i+0)+6*j+0)+=*(*(Kp+2*i+0)+2*j+0)*A*prate;
		*(*(Ke+6*i+0)+6*j+1)+=*(*(Kp+2*i+0)+2*j+1)*A*prate;
		*(*(Ke+6*i+1)+6*j+0)+=*(*(Kp+2*i+1)+2*j+0)*A*prate;
		*(*(Ke+6*i+1)+6*j+1)+=*(*(Kp+2*i+1)+2*j+1)*A*prate;
	  }
	}

	freematrix(Bpt,6);
	freematrix(DpBp,3);
	freematrix(Kp,6);
  }

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

  /*DRILL*/
  for(ii=0;ii<1;ii++)
  {
	*(*(Ke+ 5)+ 5)= 0.030*E*t*A;
	*(*(Ke+ 5)+11)=-0.015*E*t*A;
	*(*(Ke+ 5)+17)=*(*(Ke+5)+11);
	*(*(Ke+11)+ 5)=*(*(Ke+5)+11);
	*(*(Ke+11)+11)=*(*(Ke+5)+ 5);
	*(*(Ke+11)+17)=*(*(Ke+5)+11);
	*(*(Ke+17)+ 5)=*(*(Ke+5)+11);
	*(*(Ke+17)+11)=*(*(Ke+5)+11);
	*(*(Ke+17)+17)=*(*(Ke+5)+ 5);
  }

  free(w);
  freematrix(Dp,3);
  freematrix(Db,3);
  freematrix(Bp,3);
  freematrix(Bb,3);

  return Ke;
}/*assememtx*/

/*
void addkpemtx(double **Ke,struct oshell shell,double** C,double ***B)
{
  int i,j;
  double **Dp,**Bp,**Bpt,**DpBp,**Kp;
  double A = shell.sect->area;
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
  double A = shell.sect->area;
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
	t=shell.sect->t;
	A=shell.sect->area;

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






















double** consintentCshellmises（struct oshell shell, double** C, double lambda）
{
  int i,j
  double** consistentC, ** W, **H;

  return consistentC;
}

double** concistentCshellilyushin（struct oshell shell, double** C, double* lambda）
{
  int i,j
  double** consistentC, ** W, **H;

  W = (double**)malloc(6 * sizeof(double*));
  for (i = 0; i < 6; ii++)
  {
	*(W + i) = (double*)malloc(6 * sizeof(double));
  }

  ilyushin(shell, estress_try, estress, lambda, W); /*これをreturmmapに持っていくか*/

  H = matrixmatrix(W,C,6);

  consistentC = (double**)malloc(6 * sizeof(double*));
  for (i = 0; i < 6; ii++)
  {
	*(consistentC + i) = (double*)malloc(6 * sizeof(double));
	for (j = 0; j < 6; j++)
	{
		*(*(consistentC + i) + j) = *(*(H+i)+j);
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

	N1=matrixvector(H,g1,6);
	N2=matrixvector(H,g2,6);

	M11=dotproduct(g1,N1,6);
	M12=dotproduct(g1,N2,6);
	M21=M12;
	M22=dotproduct(g2,N2,6);

	det = M11*M22-M12*M21;

	for (i = 0; i < 6; ii++)
	{
	  for (j = 0; j < 6; j++)
	  {
		*(*(consistentC + i) + j) -= ( M22*(N1+i)*(N1+j) - M12*(N1+i)*(N2+j) - M21*(N2+i)*(N1+j) + M22*(N2+i)*(N2+j) )/det;
	  }
	}

	free(g1);
	free(g2);
	free(N1);
	free(N2);
  }
  else if(*(lambda+0)>0)
  {
	g1 = (double**)malloc(6 * sizeof(double*));
	*(g1+0) = a[0]/e[0] + a[3]/e[2];
	*(g1+0) = a[1]/e[0] + a[4]/e[2];
	*(g1+0) = a[2]/e[0] + a[5]/e[2];
	*(g1+0) = a[0]/e[2] + a[3]/e[1];
	*(g1+0) = a[1]/e[2] + a[4]/e[1];
	*(g1+0) = a[2]/e[2] + a[5]/e[1];

	N1=matrixvector(H,g1,6);
	M11=dotproduct(g1,N1,6);

	for (i = 0; i < 6; ii++)
	{
	  for (j = 0; j < 6; j++)
	  {
		*(*(consistentC + i) + j) -= (N1+i)*(N1+j)/M11;
	  }
	}

	free(g1);
	free(N1);
  }
  else(*(lambda+1)>0)
  {
	g2 = (double**)malloc(6 * sizeof(double*));
	*(g1+0) =  a[0]/e[0] - a[3]/e[2];
	*(g1+0) =  a[1]/e[0] - a[4]/e[2];
	*(g1+0) =  a[2]/e[0] - a[5]/e[2];
	*(g1+0) = -a[0]/e[2] + a[3]/e[1];
	*(g1+0) = -a[1]/e[2] + a[4]/e[1];
	*(g1+0) = -a[2]/e[2] + a[5]/e[1];

	N2=matrixvector(H,g2,6);
	M22=dotproduct(g2,N2,6);

	for (i = 0; i < 6; ii++)
	{
	  for (j = 0; j < 6; j++)
	  {
		*(*(consistentC + i) + j) -= (N2+i)*(N2+j)/M22;
	  }
	}

	free(g2);
	free(N2);
  }

  return consistentC;
}


double* ilyushin(struct oshell shell, double* estress_try, double* estress, double* lambda, double* W)
{
  int i,j;
  double* qstress,* ostress;
  double* f;
  double fn,fm;
  double fn2,fm2,fnfm;
  double E,t,A,I;
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
  /*[Q] = 1/√2| 1  1  0 |  */
  /*          | 0  0 √2 |  */

  /*[Q]diag[a,b,c][Q]^T  */
  /*    | a+b a-b   0 |  */
  /*=1/2| a+b a-b   0 |  */
  /*    |   0   0  2c |  */

  /*{σtry-αtry}=([I]+2Δλ[C][A]){σ-α}=[W]^-1{σ-α} */
  /*{σ-α}=[W]{σtry-αtry}                         */
  /*[W]^-1=(diag[Q,Q])[O](diag[Q^T,Q^T])         */
  /*[W]=(diag[Q,Q])[O]^-1(diag[Q^T,Q^T])         */

  E=shell.E;
  t=shell.area;

  fn=fy*t;
  fm=fy*pow(t,3)/4.0;

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




  /*(diag[Q^T,Q^T]){σtry-αtry}*/
  *(qstress+0)=( *(estress_try+0)+*(estress_try+1))/sqrt(2.0);
  *(qstress+1)=(-*(estress_try+0)+*(estress_try+1))/sqrt(2.0);
  *(qstress+2)=  *(estress_try+2);
  *(qstress+3)=( *(estress_try+3)+*(estress_try+4))/sqrt(2.0);
  *(qstress+4)=(-*(estress_try+3)+*(estress_try+4))/sqrt(2.0);
  *(qstress+5)=  *(estress_try+5);


  *(ostress+0) = Oinv[0][0][0]**(qstress+0) + Oinv[0][1][0]**(qstress+3);
  *(ostress+1) = Oinv[0][0][1]**(qstress+1) + Oinv[0][1][1]**(qstress+4);
  *(ostress+2) = Oinv[0][0][2]**(qstress+2) + Oinv[0][1][2]**(qstress+5);
  *(ostress+3) = Oinv[1][0][0]**(qstress+0) + Oinv[1][1][0]**(qstress+3);
  *(ostress+4) = Oinv[1][0][1]**(qstress+1) + Oinv[1][1][0]**(qstress+4);
  *(ostress+5) = Oinv[1][0][2]**(qstress+2) + Oinv[1][1][0]**(qstress+5);

  if(estress!=NULL)
  {
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


  /*f={σtry-αtry}^T[W]^T[A][W]{σtry-αtry}*/
  *(f+0) = p[0] * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 + *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + p[1] * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 + *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + p[2] * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 + *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm))
		 - 1.0;
  *(f+1) = p[0] * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 - *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + p[1] * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 - *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + p[2] * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 - *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm))
		 - 1.0;

  return f;
}


double* returnmapilyushin(struct oshell shell, double* estress)
{
   int i,j;
   double* lambda12,* lambda1,* lambda2,* lambdaeps;
   double* lambda;
   double* estress_try;
   double residual;
   double tolerance = 1.0e-3;
   int FLAG;

   /*FLAG = 0 ; BOTH SURFACES ARE ACTIVE*/
   /*FLAG = 1 ; SURFACE 1 IS ACTIVE     */
   /*FLAG = 2 ; SURFACE 2 IS ACTIVE     */


   lambda = (double*)malloc(2 * sizeof(double));

   /*YIELD JUDGEMENT.*/
   f=ilyushin(shell, estress, NULL, NULL);/*ILYUSHIN'S YIELD FUNCTION.{σ-α}^T[A]{σ-α}*/

   if( *(f+0)>0 || *(f+1)>0 )/*IF YIELD SURFACE ACTIVE.*/
   {
	   dlambda = (double*)malloc(2 * sizeof(double));
	   lambda12 = (double*)malloc(2 * sizeof(double));
	   lambda1  = (double*)malloc(2 * sizeof(double));
	   lambda2  = (double*)malloc(2 * sizeof(double));
	   lambdaeps = (double*)malloc(2 * sizeof(double));

	   estress_try = (double*)malloc(6 * sizeof(double));

	   /*RETURN-MAPPING PROCEDURE USING NEWTON-RAPTHON METHOD FOR ILYUSHIN'S YIELD CONDITION.*/
		for(i=0;i<2;i++)
		{
		  *(lambda12+i)=0.0;
		  *(lambda1 +i)=0.0;
		  *(lambda2 +i)=0.0;
		}
		for(i=0;i<6;i++)
		{
		  *(estress_try+i)=*(estress+i);
		}

		FLAG = 0;

		while(residual<tolerance)/*BOTH YIELD SURFACES ARE ACTIVE.*/
		{
		  residual = 0.0;
		  /*YIELD FUNCTION*/
		  f=ilyushin(shell,estress_try,NULL,lambda12);
		  /*SENSITIVITY*/
		  for(i=0;i<2;i++)
		  {
			for(j=0;j<2;j++)
			{
			  if(i==j) *(lambdaeps+i)=*(lambda12+i)+eps;
			  else     *(lambdaeps+i)=*(lambda12+i);
			}
			feps=ilyushin(shell,estress_try,NULL,lambdaeps);
			for(j=0;j<2;j++)
			{
			  *(*(K+j)+i)=(*(feps+j)-*(f+j))/eps;
			}
			free(feps);
		  }
		  free(f);

		  f11=(f11-f1)/eps;
		  f12=(f12-f1)/eps;
		  f21=(f21-f2)/eps;
		  f22=(f22-f2)/eps;

		  det=f11*f22-f12*f21;
		  *(dlambda + 0) = -( f22*f1-f12*f2)/det;
		  *(dlambda + 1) = -(-f21*f1+f11*f2)/det;

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
			while(<tolerance)/*YIELD SURFACE 1 IS ACTIVE.*//*f1=0 && f2<0 && lambda1>0 && lambda2=0*/
			{
			  /*YIELD FUNCTION*/
			  f=ilyushin(shell,estress_try,NULL,lambda1);
			  /*SENSITIVITY*/
			  *(lambdaeps+0)=*(lambda1+0)+eps;
			  *(lambdaeps+1)=0.0;

			  feps=ilyushin(shell,estress_try,NULL,lambdaeps);
			  *(feps+0)=(*(feps+0)-*(f+0))/eps;

			  *(lambda1+0) += -*(f+0) / *(feps+0);

			  free(f);
			  free(feps);
			}
		}
		if(FLAG == 2)
		{
			while(<tolerance)/*YIELD SURFACE 2 IS ACTIVE.*//*f1=0 && f2<0 && lambda1>0 && lambda2=0*/
			{
			  /*YIELD FUNCTION*/
			  f=ilyushin(shell,estress_try,NULL,lambda2);
			  /*SENSITIVITY*/
			  *(lambdaeps+1)=*(lambda2+1)+eps;
			  *(lambdaeps+0)=0.0;

			  feps=ilyushin(shell,estress_try,NULL,lambdaeps);
			  *(feps+1)=(*(feps+1)-*(f+1))/eps;

			  *(lambda2+1) += -*(f+1) / *(feps+1);

			  free(f);
			  free(feps);
			}
		}

		/*CHOOSE YIELD CONDITION & UPDATE LAMBDA*/
		for(i=0;i<2;i++)
		{
		  if()
		  {
			*(lambda+i)=*(lambda12+i);
		  }
		  else if()
		  {
			*(lambda+i)=*(lambda1+i);
		  }
		  else
		  {
			*(lambda+i)=*(lambda2+i);
		  }
		}

		/*UPDATE STRESS*/

		//f=ilyushin(estress_try,estress,lambda);

		free(estress_try);
		free(lambdaeps);
		free(lambda12);
		free(lambda1);
		free(lambda2);
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
