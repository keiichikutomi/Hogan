/*FUNCTIONS FOR PLASTICITY*/




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








		/*[P]=2[A]*/
		/*{dε_p}=dλ*{df/dσ}=dλ[P]{σ}*/













double** Cmises（struct oshell shell, double lambda）
{

}



double** concistentCilyushin（struct oshell shell, double** C, double* lambda）
{
  int i,j
  double** consistentC, ** W;




  W = (double**)malloc(6 * sizeof(double*));
  for (i = 0; i < 6; ii++)
  {
	*(W + i) = (double*)malloc(6 * sizeof(double));
  }


  ilyushin(shell, estress_try, estress, lambda, W); /*これをreturmmapに持っていくか*/

  H = matrixmatrix(W,C,6);

	A1 = (double**)malloc(6 * sizeof(double*));
	for (ii = 0; ii < 6; ii++)
	{
		*(A1 + ii) = (double*)malloc(6 * sizeof(double));
		for (jj = 0; jj < 6; jj++)
		{
			*(*(A + ii) + jj) = 0.0;
		}
	}

	*(*(A1+0)+0)=1.0/fn2;
	*(*(A1+0)+1)=-0.5/fn2;
	*(*(A1+1)+0)=*(*(A1+0)+1);
	*(*(A1+1)+1)=*(*(A1+0)+0);
	*(*(A1+2)+2)=3.0/fn2;

	*(*(A1+0)+3)=1.0/(2.0*sqrt(3.0)*fnfm);
	*(*(A1+0)+4)=-0.5/(2.0*sqrt(3.0)*fnfm);
	*(*(A1+1)+3)=*(*(A1+0)+4);
	*(*(A1+1)+4)=*(*(A1+0)+3);
	*(*(A1+2)+5)=3.0/(2.0*sqrt(3.0)*fnfm);

	*(*(A1+3)+0)=*(*(A1+0)+3);
	*(*(A1+3)+1)=*(*(A1+1)+3);
	*(*(A1+4)+0)=*(*(A1+0)+4);
	*(*(A1+4)+1)=*(*(A1+1)+4);
	*(*(A1+5)+2)=3.0/(2.0*sqrt(3.0)*fnfm);

	*(*(A1+3)+3)=1.0/fm2;
	*(*(A1+3)+4)=-0.5/fm2;
	*(*(A1+4)+3)=*(*(A1+3)+4);
	*(*(A1+4)+4)=*(*(A1+3)+3);
	*(*(A1+5)+5)=3.0/fm2;



	A2 = (double**)malloc(6 * sizeof(double*));
	for (ii = 0; ii < 6; ii++)
	{
		*(A2 + ii) = (double*)malloc(6 * sizeof(double));
		for (jj = 0; jj < 6; jj++)
		{
			if((ii<3 && jj<3) || (ii>2 && jj>2)) *(*(A2 + ii) + jj) = *(*(A1 + ii) + jj);
			else *(*(A2 + ii) + jj) = -*(*(A1 + ii) + jj);
		}
	}

  g=matrixvector(A,estress);

  Hg=matrixvector(H,g);
  gHg=dotproduct(g,Hg,6);

  consistentC = (double**)malloc(6 * sizeof(double*));
  for (i = 0; i < 6; ii++)
  {
	*(consistentC + i) = (double*)malloc(6 * sizeof(double));
	for (j = 0; j < 6; j++)
	{
	  *(*(consistentC + i) + j) = *(*(H+i)+j) + *(Hg+i)**(Hg+j)/(beta+gHg);
	}
  }








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

   lambda = (double*)malloc(2 * sizeof(double));

   /*YIELD JUDGEMENT.*/
   f=ilyushin(shell, estress, NULL, NULL);/*ILYUSHIN'S YIELD FUNCTION.{σ-α}^T[A]{σ-α}*/

   if( *(f+0)>0 || *(f+1)>0 )/*IF YIELD SURFACE ACTIVE.*/
   {

	   lambda12 = (double*)malloc(2 * sizeof
	   (double));
	   lambda1 = (double*)malloc(2 * sizeof(double));
	   lambda2 = (double*)malloc(2 * sizeof(double));
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

		while(<tolerance)/*BOTH YIELD SURFACES ARE ACTIVE.*/
		{
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
		  lambda1 += -( f22*f1-f12*f2)/det;
		  lambda2 += -(-f21*f1+f11*f2)/det;
		}

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
