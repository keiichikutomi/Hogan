


/*TYPE |      JOINT        : DESCRIPTION                                                  */
/*  1  | REVOLUTE JOINT    : 1 DIRECTION ROTATION, NO TRANSITION.          1 VECTOR INPUT.*/
/*  2  | SPHERICAL JOINT   : FREE ROTATION, NO TRANSITION.      		   0 VECTOR INPUT.*/
/*  3  | PRISMATIC JOINT   : 1 DIRECTION TRANSITION, NO ROTATION.		   1 VECTOR INPUT.*/
/*  4  | CYLINDRICAL JOINT : 1 DIRECTION ROTATION, 1 DIRECTION TRANSITION. 2 VECTOR INPUT.*/
/*  5  | UNIVERSAL JOINT   : 2 DIRECTION ROTATION, NO TRANSITION.          2 VECTOR INPUT.*/

/*
for (i = 0; i < nconstraint; i++)
{
	msize += (constraints+i)->neq;
}
*/

void assemconstraint(struct oconstraint* constraints, int nconstraint, long int* constraintmain,
					 struct gcomponent* gmtx, double* gvct,
					 double* iform, double* ddisp, int msize)
{
	char str[500];
	int i,j,ii,jj;
	int* loffset;

	int nnod,type,neq,leq;

	double dot;
	double *cross;

	double* initaxis;
	double* axis1,* axis2;
	double* rvct1,* rvct2;
	double** rmtx1,** rmtx2;

	double* cvct;
	double** cmtx;

	for (i = 1; i <= nconstraint; i++)
	{
		neq = (constraints+i-1)->neq;
		leq = (constraints+i-1)->leq;
		type = (constraints+i-1)->type;


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




		if(type == 1)
		{
			nnod = (constraints+i-1)->nnod;

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

			rvct1 = (double*)malloc(3 * sizeof(double));
			rvct2 = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(rvct1 + ii) = *(ddisp + *(loffset + 3 + ii));
				*(rvct2 + ii) = *(ddisp + *(loffset + 9 + ii));
			}
			rmtx1 = rotationmtx(rvct1);
			rmtx2 = rotationmtx(rvct2);

			initaxis = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(initaxis + ii) = (constraints+i-1)->axis[ii];
			}

			axis1 = matrixvector(rmtx1,initaxis,3);
			axis2 = matrixvector(rmtx2,initaxis,3);



			dot = dotproduct(axis1,axis2,3);
			*(cvct+0) = *(ddisp + *(loffset + 0)) - *(ddisp + *(loffset + 6));
			*(cvct+1) = *(ddisp + *(loffset + 1)) - *(ddisp + *(loffset + 7));
			*(cvct+2) = *(ddisp + *(loffset + 2)) - *(ddisp + *(loffset + 8));
			*(cvct+3) = dot - 1.0;

			cross = crossproduct(axis1,axis2);
			*(*(cmtx+0)+0) = 1.0; *(*(cmtx+0)+6) = -1.0;
			*(*(cmtx+1)+1) = 1.0; *(*(cmtx+1)+7) = -1.0;
			*(*(cmtx+2)+2) = 1.0; *(*(cmtx+2)+8) = -1.0;
			for (ii = 0; ii < 3; ii++)
			{
				*(*(cmtx+3)+3+ii) =  *(cross+ii);
				*(*(cmtx+3)+9+ii) = -*(cross+ii);
			}

			for (ii = 0; ii < neq; ii++)
			{
				*(gvct + msize + leq + ii) = *(cvct + ii);
			}


			for (ii = 0; ii < neq; ii++)
			{
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					if(*(*(cmtx+ii)+jj) != 0)
					{
						gwrite(gmtx , msize + leq + ii + 1, *(constraintmain + *(loffset + jj)) + 1, *(*(cmtx+ii)+jj));
					}
				}
			}

			free(loffset);
			free(rvct1);
			free(rvct2);
			freematrix(rmtx1,3);
			freematrix(rmtx2,3);
			free(axis1);
			free(axis2);
			free(cross);
			free(cvct);
			freematrix(cmtx,neq);
		}

	}
	return;
}




