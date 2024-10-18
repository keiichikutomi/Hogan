/*FUNCTIONS FOR PLASTICITY*/









C = (double**)malloc(6 * sizeof(double*));
for (ii = 0; ii < 6; ii++)
{
	*(C + ii) = (double*)malloc(6 * sizeof(double));
	for (jj = 0; jj < 6; jj++)
	{
		*(*(C + ii) + jj) = 0.0;
	}
}

*(*(C+0)+0)=E*t/(1-poi*poi);
*(*(C+0)+1)=poi**(*(C+0)+0));
*(*(C+0)+2)=0.0;
*(*(C+1)+0)=*(*(C+0)+1);
*(*(C+1)+1)=*(*(C+0)+0);
*(*(C+1)+2)=0.0;
*(*(C+2)+0)=0.0;
*(*(C+2)+1)=0.0;
*(*(C+2)+2)=0.5*E*t/(1+poi);

*(*(C+3)+3)=E*pow(t,3)/(12.0*(1-poi*poi));
*(*(C+3)+4)=poi**(*(C+0)+0));
*(*(C+3)+5)=0.0;
*(*(C+4)+3)=*(*(C+3)+4);
*(*(C+4)+4)=*(*(C+3)+3);
*(*(C+4)+5)=0.0;
*(*(C+5)+3)=0.0;
*(*(C+5)+4)=0.0;
*(*(C+5)+5)=0.5*E*pow(t,3)/(12.0*(1+poi));


		/*[H]^-1=[C]^-1+dƒÉ*[P]*/
		Cinv = (double**)malloc(6 * sizeof(double*));
		for (ii = 0; ii < 6; ii++)
		{
			*(Cinv + ii) = (double*)malloc(6 * sizeof(double));
			for (jj = 0; jj < 6; jj++)
			{
				*(*(Cinv + ii) + jj) = 0.0;
			}
		}

		*(*(H+0)+0)=1.0/(E*t) 					+ 2.0*lambda1**(*(A1+0)+0);
		*(*(H+0)+1)=-poi/(E*t) 					+ 2.0*lambda1**(*(A1+0)+1);
		*(*(H+1)+0)=*(*(H+0)+1);
		*(*(H+1)+1)=*(*(H+0)+0);
		*(*(H+2)+2)=2.0*(1.0+poi)/(E*t) 		    + 2.0*lambda1**(*(A1+2)+2);

		*(*(H+3)+3)=12.0/(E*pow(t,3))			    + 2.0*lambda1**(*(A1+3)+3);
		*(*(H+3)+4)=-12.0*poi/(E*pow(t,3)) 		+ 2.0*lambda1**(*(A1+3)+4);
		*(*(H+4)+3)=*(*(H+3)+4);
		*(*(H+4)+4)=*(*(H+3)+3);
		*(*(H+5)+5)=24.0*(1.0+poi)/(E*pow(t,3))   + 2.0*lambda1**(*(A1+5)+5);

		/*[P]=2[A]*/
		/*{dƒÃ_p}=dƒÉ*{df/dƒÐ}=dƒÉ[P]{ƒÐ}*/






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




		  Hg=matrixvector(H,g);
		  gHg=dotproduct(g,Hg,6);



		  consistentC = (double**)malloc(6 * sizeof(double*));
		  for (ii = 0; ii < 6; ii++)
		  {
			*(consistentC + ii) = (double*)malloc(6 * sizeof(double));
			for (jj = 0; jj < 6; jj++)
			{
			  *(*(consistentC + ii) + jj) = *(*(H+ii)+jj) + *(Hg+ii)**(Hg+jj)/(beta+gHg);
			}
		  }



double* ilyushin(struct oshell shell, double* estress_try, double* estress, double* lambda)
{
  double* qstress,* ostress;
  double f;
  double E,t;
  double lambda1,lambda2;

  double c1,c2,c3;
  double O00,O11,O22,O33,O44,O55,O03,O14,O25,O30,O41,O52;
  double Oinv00,Oinv11,Oinv22,Oinv33,Oinv44,Oinv55,Oinv03,Oinv14,Oinv25,Oinv30,Oinv41,Oinv52;
  double det03,det14,det25;


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
  E=shell.E;
  t=shell.area;


  /*           | 1 -1  0 |  */
  /* [Q] = 1/ã2| 1  1  0 |  */
  /*           | 0  0 ã2 |  */

  /*(diag[Q^T,Q^T]){ƒÐtry-ƒ¿try}*/
  *(qstress+0)=( *(estress_try+0)+*(estress_try+1))/sqrt(2.0);
  *(qstress+1)=(-*(estress_try+0)+*(estress_try+1))/sqrt(2.0);
  *(qstress+2)=  *(estress_try+2);
  *(qstress+3)=( *(estress_try+3)+*(estress_try+4))/sqrt(2.0);
  *(qstress+4)=(-*(estress_try+3)+*(estress_try+4))/sqrt(2.0);
  *(qstress+5)=  *(estress_try+5);

  /*[W]=(diag[Q,Q])[O]^-1(diag[Q^T,Q^T])*/

  /*[O]*/



  /*{ƒÐtry-ƒ¿try}=([I]+2ƒ¢ƒÉ[C][A]){ƒÐ-ƒ¿}              */
  /*           =(diag[Q,Q])[O](diag[Q^T,Q^T]){ƒÐ-ƒ¿}*/

  /*([C]^-1+2ƒ¢ƒÉ[A])*/
  fn=fy*t;
  fm=fy*pow(t,3)/4.0;

  fn2=pow(fn,2);
  fm2=pow(fm,2);
  fnfm=fn*fm;

  c1 = (lambda1 + lambda2) / fn2;
  c2 = (lambda1 + lambda2) / fm2;
  c3 = (lambda1 - lambda2) /(2.0*sqrt(3.0)*fnfm);

  O00 = E*t/(1-v)               /*+ 2.0/3.0*Hkin*/;
  O11 = 3.0*E*t/(1+v)           /*+ 2.0*Hkin*/    ;
  O22 = 3.0*E*t/(1+v)           /*+ 4.0*Hkin*/    ;
  O33 = E*pow(t,3)/(12.0*(1-v)) /*+ 2.0/3.0*Hkin*/;
  O44 = E*pow(t,3)/( 4.0*(1+v)) /*+ 2.0*Hkin*/    ;
  O55 = E*pow(t,3)/( 4.0*(1+v)) /*+ 4.0*Hkin*/    ;
  O03 = O00;
  O14 = O11;
  O25 = O22;
  O30 = O33;
  O41 = O44;
  O52 = O55;

  O00 = 1.0 + c1*O00;
  O11 = 1.0 + c1*O11;
  O22 = 1.0 + c1*O22;
  O33 = 1.0 + c2*O33;
  O44 = 1.0 + c2*O44;
  O55 = 1.0 + c2*O55;
  O03 *= c3;
  O14 *= c3;
  O25 *= c3;
  O30 *= c3;
  O41 *= c3;
  O52 *= c3;


  /*Oinv:[O]^-1*/
  det03 = O00*O33-O03*O30;
  det14 = O11*O44-O14*O41;
  det25 = O22*O55-O25*O52;

  Oinv00 =  O33 / det03;
  Oinv11 =  O44 / det14;
  Oinv22 =  O55 / det25;
  Oinv33 =  O00 / det03;
  Oinv44 =  O11 / det14;
  Oinv55 =  O22 / det25;
  Oinv03 = -O03 / det03;
  Oinv14 = -O14 / det14;
  Oinv25 = -O25 / det25;
  Oinv30 = -O30 / det03;
  Oinv41 = -O41 / det14;
  Oinv52 = -O52 / det25;


  *(ostress+0)= Oinv00**(qstress+0) + Oinv03**(qstress+3);
  *(ostress+1)= Oinv11**(qstress+1) + Oinv14**(qstress+4);
  *(ostress+2)= Oinv22**(qstress+2) + Oinv25**(qstress+5);
  *(ostress+3)= Oinv30**(qstress+0) + Oinv33**(qstress+3);
  *(ostress+4)= Oinv41**(qstress+1) + Oinv44**(qstress+4);
  *(ostress+5)= Oinv52**(qstress+2) + Oinv55**(qstress+5);

  if(estress!=NULL)
  {
	*(estress+0)=(*(ostress+0)-*(ostress+1))/sqrt(2.0);
	*(estress+1)=(*(ostress+0)+*(ostress+1))/sqrt(2.0);
	*(estress+2)= *(ostress+2);
	*(estress+3)=(*(ostress+3)-*(ostress+4))/sqrt(2.0);
	*(estress+4)=(*(ostress+3)+*(ostress+4))/sqrt(2.0);
	*(estress+5)= *(ostress+5);
  }



  /*f={ƒÐtry-ƒ¿try}^T[W]^T[A][W]{ƒÐtry-ƒ¿try}*/

  *(f+0) = 0.5 * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 + *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + 1.5 * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 + *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + 3.0 * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 + *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm))
		 - 1.0;
  *(f+1) = 0.5 * (*(ostress+0)**(ostress+0)/fn2 + *(ostress+3)**(ostress+3)/fm2 - *(ostress+0)**(ostress+3)/(sqrt(3.0)*fnfm))
		 + 1.5 * (*(ostress+1)**(ostress+1)/fn2 + *(ostress+4)**(ostress+4)/fm2 - *(ostress+1)**(ostress+4)/(sqrt(3.0)*fnfm))
		 + 3.0 * (*(ostress+2)**(ostress+2)/fn2 + *(ostress+5)**(ostress+5)/fm2 - *(ostress+2)**(ostress+5)/(sqrt(3.0)*fnfm))
		 - 1.0;

  return f;
}


double* returnmapIlyushin(struct oshell shell, double* estress)
{
   int i,j;
   double* lambda12,* lambda1,* lambda2,* lambdaeps;
   double* lambda;
   double* estress_try;

   lambda = (double*)malloc(2 * sizeof(double));

   /*YIELD JUDGEMENT.*/
   f=ilyushin(shell, estress, NULL, NULL);/*ILYUSHIN'S YIELD FUNCTION.{ƒÐ-ƒ¿}^T[A]{ƒÐ-ƒ¿}*/

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
		  f=ilyushin(estress_try,NULL,lambda12);
		  /*SENSITIVITY*/
		  for(i=0;i<2;i++)
		  {
			for(j=0;j<2;j++)
			{
			  if(i==j) *(lambdaeps+i)=*(lambda12+i)+eps;
			  else     *(lambdaeps+i)=*(lambda12+i);
			}
			feps=ilyushin(estress_try,NULL,lambdaeps);
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
		  f=ilyushin(estress_try,NULL,lambda1);
		  /*SENSITIVITY*/
		  *(lambdaeps+0)=*(lambda1+0)+eps;
		  *(lambdaeps+1)=0.0;

		  feps=ilyushin(estress_try,NULL,lambdaeps);
		  *(feps+0)=(*(feps+0)-*(f+0))/eps;

		  *(lambda1+0) += -*(f+0) / *(feps+0);

		  free(f);
		  free(feps);
		}

		while(<tolerance)/*YIELD SURFACE 2 IS ACTIVE.*//*f1=0 && f2<0 && lambda1>0 && lambda2=0*/
		{
		  /*YIELD FUNCTION*/
		  f=ilyushin(estress_try,NULL,lambda2);
		  /*SENSITIVITY*/
		  *(lambdaeps+1)=*(lambda2+1)+eps;
		  *(lambdaeps+0)=0.0;

		  feps=ilyushin(estress_try,NULL,lambdaeps);
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

		f=ilyushin(estress_try,estress,lambda);

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
