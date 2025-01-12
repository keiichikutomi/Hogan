


/*TYPE |      JOINT        : DESCRIPTION                                                  */
/*  1  | REVOLUTE JOINT    : 1 DIRECTION ROTATION, NO TRANSITION.          1 VECTOR INPUT.*/
/*  2  | SPHERICAL JOINT   : FREE ROTATION, NO TRANSITION.      		   0 VECTOR INPUT.*/
/*  3  | PRISMATIC JOINT   : 1 DIRECTION TRANSITION, NO ROTATION.		   1 VECTOR INPUT.*/
/*  4  | CYLINDRICAL JOINT : 1 DIRECTION ROTATION, 1 DIRECTION TRANSITION. 2 VECTOR INPUT.*/
/*  5  | UNIVERSAL JOINT   : 2 DIRECTION ROTATION, NO TRANSITION.          2 VECTOR INPUT.*/



void assemconstraint(struct oconstraint* constraints, int nconstraint, long int* constraintmain,
					 struct gcomponent* gmtx,
					 double* iform, double* ddisp, double* lambda, int msize)
{
	char str[500];
	int i,j,ii,jj;
	int* loffset;

	int nnod,type,neq,leq;

	double dot1 = 0.0;
	double dot2 = 0.0;
	double *cross1,*cross2;

	double* initaxis1,* initaxis2,* initaxis3;
	double* axis1,* axis2,* axis3;
	double* rvct1,* rvct2;
	double** rmtx1,** rmtx2;

	double* cvct;
	double** cmtx;

	for (i = 1; i <= nconstraint; i++)
	{
		neq = (constraints+i-1)->neq;
		leq = (constraints+i-1)->leq;
		type = (constraints+i-1)->type;
		nnod = (constraints+i-1)->nnod;


			  /*ASSEMBLAGE CONSTRAINT MATRIX.*/
	  /*
	  for(i=1;i<=nconstraint;i++)
	  {
		for (ii=1;ii<=msize;ii++)
		{
		  if(*(*(constraintvec+i-1)+ii-1) != 0 || *(*(constraintvec+i-1)+ii-1) != NULL)
		  {
			gdata=*(*(constraintvec+i-1)+ii-1);
			gwrite(gmtx,msize+i,ii,gdata);
			comps++;
		  }
		}
		*(gvct+msize+i-1) = *(constraintval+i-1);
		(confs+msize+i-1)->iconf = (signed char)0;
		(confs+msize+i-1)->value = *(constraintval+i-1);
	  }
	  */




		if(type == 1)/*REVOLUTE*/
		{
			cmtx = (double**)malloc(neq * sizeof(double*));
			for(ii = 0 ; ii < neq ; ii++)
			{
				*(cmtx+ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6*nnod; jj++)
				{
					*(*(cmtx+ii)+jj) = 0.0;
				}
			}

			loffset = (int*)malloc(6 * nnod * sizeof(int));
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(loffset + (6 * ii + jj)) = 6 * ((constraints+i-1)->node[ii]->loff) + jj;
				}
			}

			rvct1 = (double*)malloc(3 * sizeof(double));
			rvct2 = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(rvct1 + ii) = *(ddisp + *(loffset + 3 + ii));
				*(rvct2 + ii) = *(ddisp + *(loffset + 9 + ii));
			}
			rmtx1 = rotationmtx(rvct1);
			rmtx2 = rotationmtx(rvct2);

			initaxis1 = (double*)malloc(3 * sizeof(double));
			initaxis2 = (double*)malloc(3 * sizeof(double));
			initaxis3 = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(initaxis1 + ii) = (constraints+i-1)->axis[0][ii];
				*(initaxis2 + ii) = (constraints+i-1)->axis[1][ii];
				*(initaxis3 + ii) = (constraints+i-1)->axis[2][ii];
			}

			axis1 = matrixvector(rmtx1,initaxis1,3);
			axis2 = matrixvector(rmtx1,initaxis2,3);
			axis3 = matrixvector(rmtx2,initaxis3,3);

			cross1 = crossproduct(axis1,axis3);
			cross2 = crossproduct(axis2,axis3);


			*(*(cmtx+0)+0) = 1.0; *(*(cmtx+0)+6) = -1.0;
			*(*(cmtx+1)+1) = 1.0; *(*(cmtx+1)+7) = -1.0;
			*(*(cmtx+2)+2) = 1.0; *(*(cmtx+2)+8) = -1.0;
			for (ii = 0; ii < 3; ii++)
			{
				*(*(cmtx+3)+3+ii) =  *(cross1+ii);
				*(*(cmtx+3)+9+ii) = -*(cross1+ii);
				*(*(cmtx+4)+3+ii) =  *(cross2+ii);
				*(*(cmtx+4)+9+ii) = -*(cross2+ii);
			}


			for (ii = 0; ii < neq; ii++)
			{
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					if(*(*(cmtx+ii)+jj) != 0)
					{
						gwrite(gmtx , msize+leq+ii+1, *(constraintmain + *(loffset + jj)) + 1, *(*(cmtx+ii)+jj));
					}
				}
			}

			free(loffset);
			free(rvct1);
			free(rvct2);
			freematrix(rmtx1,3);
			freematrix(rmtx2,3);
			free(initaxis1);
			free(initaxis2);
			free(initaxis3);
			free(axis1);
			free(axis2);
			free(axis3);
			free(cross1);
			free(cross2);
			freematrix(cmtx,neq);
		}
		if(type == 2)/*SPHERICAL*/
		{
			cmtx = (double**)malloc(neq * sizeof(double*));
			for(ii = 0 ; ii < neq ; ii++)
			{
				*(cmtx+ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6*nnod; jj++)
				{
					*(*(cmtx+ii)+jj) = 0.0;
				}
			}

			loffset = (int*)malloc(6 * nnod * sizeof(int));
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(loffset + (6 * ii + jj)) = 6 * ((constraints+i-1)->node[ii]->loff) + jj;
				}
			}

			*(*(cmtx+0)+0) = 1.0; *(*(cmtx+0)+6) = -1.0;
			*(*(cmtx+1)+1) = 1.0; *(*(cmtx+1)+7) = -1.0;
			*(*(cmtx+2)+2) = 1.0; *(*(cmtx+2)+8) = -1.0;


			for (ii = 0; ii < neq; ii++)
			{
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					if(*(*(cmtx+ii)+jj) != 0)
					{
						gwrite(gmtx , msize+leq+ii+1, *(constraintmain + *(loffset + jj)) + 1, *(*(cmtx+ii)+jj));
					}
				}
			}

			free(loffset);

			freematrix(cmtx,neq);
		}


	}
	return;
}




void constraintstress(struct oconstraint* constraints, int nconstraint, long int* constraintmain,
					  double* iform, double* ddisp, double* lambda, double* finternal, double* constraintvct)
{
	char str[500];
	int i,j,ii,jj;
	int* loffset;

	int nnod,type,neq,leq;

	double dot1 = 0.0;
	double dot2 = 0.0;
	double *cross1,*cross2;

	double* initaxis1,* initaxis2,* initaxis3;
	double* axis1,* axis2,* axis3;
	double* rvct1,* rvct2;
	double** rmtx1,** rmtx2;
	double* ginternal;

	double* cvct;
	double** cmtx;

	for (i = 1; i <= nconstraint; i++)/*REVOLUTE*/
	{

		neq = (constraints+i-1)->neq;
		leq = (constraints+i-1)->leq;
		type = (constraints+i-1)->type;
		nnod = (constraints+i-1)->nnod;

		if(type == 1)
		{
			cvct = (double*)malloc(neq * sizeof(double));
			cmtx = (double**)malloc(neq * sizeof(double*));
			for(ii = 0 ; ii < neq ; ii++)
			{
				*(cmtx+ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6*nnod; jj++)
				{
					*(*(cmtx+ii)+jj) = 0.0;
				}
			}


			loffset = (int*)malloc(6 * nnod * sizeof(int));
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(loffset + (6 * ii + jj)) = 6 * ((constraints+i-1)->node[ii]->loff) + jj;
				}
			}

			//sprintf(str,"%d %d %d %d %d %d\n", *(loffset + 0),*(loffset + 1),*(loffset + 2),*(loffset + 6),*(loffset + 7),*(loffset + 8));
			//errormessage(str);

			rvct1 = (double*)malloc(3 * sizeof(double));
			rvct2 = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(rvct1 + ii) = *(ddisp + *(loffset + 3 + ii));
				*(rvct2 + ii) = *(ddisp + *(loffset + 9 + ii));
			}
			rmtx1 = rotationmtx(rvct1);
			rmtx2 = rotationmtx(rvct2);

			initaxis1 = (double*)malloc(3 * sizeof(double));
			initaxis2 = (double*)malloc(3 * sizeof(double));
			initaxis3 = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(initaxis1 + ii) = (constraints+i-1)->axis[0][ii];
				*(initaxis2 + ii) = (constraints+i-1)->axis[1][ii];
				*(initaxis3 + ii) = (constraints+i-1)->axis[2][ii];
			}

			axis1 = matrixvector(rmtx1,initaxis1,3);
			axis2 = matrixvector(rmtx1,initaxis2,3);
			axis3 = matrixvector(rmtx2,initaxis3,3);


			dot1 = dotproduct(axis1,axis3,3);
			dot2 = dotproduct(axis2,axis3,3);
			*(cvct+0) = *(ddisp + *(loffset + 0)) - *(ddisp + *(loffset + 6));
			*(cvct+1) = *(ddisp + *(loffset + 1)) - *(ddisp + *(loffset + 7));
			*(cvct+2) = *(ddisp + *(loffset + 2)) - *(ddisp + *(loffset + 8));
			*(cvct+3) = dot1;
			*(cvct+4) = dot2;

			cross1 = crossproduct(axis1,axis3);
			cross2 = crossproduct(axis2,axis3);
			*(*(cmtx+0)+0) = 1.0; *(*(cmtx+0)+6) = -1.0;
			*(*(cmtx+1)+1) = 1.0; *(*(cmtx+1)+7) = -1.0;
			*(*(cmtx+2)+2) = 1.0; *(*(cmtx+2)+8) = -1.0;
			for (ii = 0; ii < 3; ii++)
			{
				*(*(cmtx+3)+3+ii) =  *(cross1+ii);
				*(*(cmtx+3)+9+ii) = -*(cross1+ii);
				*(*(cmtx+4)+3+ii) =  *(cross2+ii);
				*(*(cmtx+4)+9+ii) = -*(cross2+ii);
			}


			ginternal = (double*)malloc(6*nnod * sizeof(double));
			for (ii = 0; ii < 6*nnod; ii++)
			{
				*(ginternal + ii) = 0.0;
				for (jj = 0; jj < neq; jj++)
				{
					*(ginternal + ii) += *(lambda+leq+jj)**(*(cmtx+jj)+ii);
				}
				*(finternal + *(constraintmain + *(loffset + ii))) += *(ginternal + ii);
			}

			for (ii = 0; ii < neq; ii++)
			{
				*(constraintvct+leq+ii) -= *(cvct+ii);
			}


			free(loffset);
			free(rvct1);
			free(rvct2);
			freematrix(rmtx1,3);
			freematrix(rmtx2,3);
			free(initaxis1);
			free(initaxis2);
			free(initaxis3);
			free(axis1);
			free(axis2);
			free(axis3);
			free(cross1);
			free(cross2);
			freematrix(cmtx,neq);
			free(cvct);
			free(ginternal);
		}
		if(type == 2)/*SPHERICAL*/
		{
			cvct = (double*)malloc(neq * sizeof(double));
			cmtx = (double**)malloc(neq * sizeof(double*));
			for(ii = 0 ; ii < neq ; ii++)
			{
				*(cmtx+ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6*nnod; jj++)
				{
					*(*(cmtx+ii)+jj) = 0.0;
				}
			}


			loffset = (int*)malloc(6 * nnod * sizeof(int));
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(loffset + (6 * ii + jj)) = 6 * ((constraints+i-1)->node[ii]->loff) + jj;
				}
			}

			*(cvct+0) = *(ddisp + *(loffset + 0)) - *(ddisp + *(loffset + 6));
			*(cvct+1) = *(ddisp + *(loffset + 1)) - *(ddisp + *(loffset + 7));
			*(cvct+2) = *(ddisp + *(loffset + 2)) - *(ddisp + *(loffset + 8));


			*(*(cmtx+0)+0) = 1.0; *(*(cmtx+0)+6) = -1.0;
			*(*(cmtx+1)+1) = 1.0; *(*(cmtx+1)+7) = -1.0;
			*(*(cmtx+2)+2) = 1.0; *(*(cmtx+2)+8) = -1.0;


			ginternal = (double*)malloc(6*nnod * sizeof(double));
			for (ii = 0; ii < 6*nnod; ii++)
			{
				*(ginternal + ii) = 0.0;
				for (jj = 0; jj < neq; jj++)
				{
					*(ginternal + ii) += *(lambda+leq+jj)**(*(cmtx+jj)+ii);
				}
				*(finternal + *(constraintmain + *(loffset + ii))) += *(ginternal + ii);
			}

			for (ii = 0; ii < neq; ii++)
			{
				*(constraintvct+leq+ii) -= *(cvct+ii);
			}


			free(loffset);
			freematrix(cmtx,neq);
			free(cvct);
			free(ginternal);
		}

	}
	return;
}


