


/*TYPE |      JOINT        : DESCRIPTION                                                  */
/*  1  | REVOLUTE JOINT    : 1 DIRECTION ROTATION, NO TRANSITION.          1 VECTOR INPUT.*/
/*  2  | SPHERICAL JOINT   : FREE ROTATION, NO TRANSITION.      		   0 VECTOR INPUT.*/
/*  3  | PRISMATIC JOINT   : 1 DIRECTION TRANSITION, NO ROTATION.		   1 VECTOR INPUT.*/
/*  4  | CYLINDRICAL JOINT : 1 DIRECTION ROTATION, 1 DIRECTION TRANSITION. 2 VECTOR INPUT.*/
/*  5  | UNIVERSAL JOINT   : 2 DIRECTION ROTATION, NO TRANSITION.          2 VECTOR INPUT.*/
/*  6  | PERIODIC BOUNDARY :                                                              */



double* extractconstraintdisplacement(struct oconstraint constraint, double* ddisp)
/*EXTRACT ELEMENT DEFORMATION{dU} FROM GLOBAL VECTOR.*/
{
	long int i, loffset;
	int n, nnod;
	double* d;

	nnod = constraint.nnod;
	d = (double*)malloc(6 * nnod * sizeof(double));
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 6; i++)
		{
			loffset = 6 * (constraint.node[n]->loff) + i;
			*(d + 6 * n + i) = *(ddisp + loffset);
		}
	}
	return d;
}/*extractconstraintdisplacement*/

void assemconstraint(struct oconstraint* constraints, int nconstraint, long int* constraintmain,
					 struct gcomponent* gmtx,
					 double* iform, double* ddisp, double* lambda, int msize)
{
	char str[500];
	int i,j,ii,jj;
	int* loffset;
	int nnod,type,neq,leq;

	double* cvct;
	double** cmtx;
	double** cmtxH;

	double gdata;

	double* initaxis1,* initaxis2,* initaxis3;
	double* axis1,* axis2,* axis3;
	double *cross1,*cross2;
	double dot1 = 0.0;
	double dot2 = 0.0;

	double *gforminit, *gform, *gdisp;

	double *rvct_m, *rvct_s, *rvct_c, *rvct_cmst;
	double **R_m, **R_s, **R_c;
	double **R_st, **R_cm, **R_cmst;
	double **H, **HR_cmst, **HR_c;


	for (i = 1; i <= nconstraint; i++)
	{
		neq = (constraints+i-1)->neq;
		leq = (constraints+i-1)->leq;
		type = (constraints+i-1)->type;
		nnod = (constraints+i-1)->nnod;


		if(type == 1)/*REVOLUTE*/
		{
			/*Jacobian*/
			cmtx = (double**)malloc(neq * sizeof(double*));
			for(ii = 0 ; ii < neq ; ii++)
			{
				*(cmtx+ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6*nnod; jj++)
				{
					*(*(cmtx+ii)+jj) = 0.0;
				}
			}

			/*Hessian*/
			cmtxH = (double**)malloc(6 * nnod * sizeof(double*));
			for(ii = 0 ; ii < 6 * nnod ; ii++)/*HESSIAN OF CONSTRAINT*/
			{
				*(cmtxH+ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					*(*(cmtxH+ii)+jj) = 0.0;
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

			gforminit = extractconstraintdisplacement(*(constraints+i-1), iform);
			gform = extractconstraintdisplacement(*(constraints+i-1), ddisp);
			gdisp = extractdeformation(gforminit, gform, nnod);


			/*Rotation Vector*/
			rvct_m = (double*)malloc(3 * sizeof(double));
			rvct_s = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(rvct_m + ii) = *(gdisp + *(loffset + 3 + ii));
				*(rvct_s + ii) = *(gdisp + *(loffset + 9 + ii));
			}
			R_m = rotationmtx(rvct_m);
			R_s = rotationmtx(rvct_s);

			initaxis1 = (double*)malloc(3 * sizeof(double));
			initaxis2 = (double*)malloc(3 * sizeof(double));
			initaxis3 = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(initaxis1 + ii) = (constraints+i-1)->axis[0][ii];
				*(initaxis2 + ii) = (constraints+i-1)->axis[1][ii];
				*(initaxis3 + ii) = (constraints+i-1)->axis[2][ii];
			}

			axis1 = matrixvector(R_m,initaxis1,3);
			axis2 = matrixvector(R_m,initaxis2,3);
			axis3 = matrixvector(R_s,initaxis3,3);

			cross1 = crossproduct(axis1,axis3);
			cross2 = crossproduct(axis2,axis3);


			for (ii = 0; ii < 3; ii++)
			{
				*(*(cmtx+0)+3+ii) =  *(cross1+ii);
				*(*(cmtx+0)+9+ii) = -*(cross1+ii);
				*(*(cmtx+1)+3+ii) =  *(cross2+ii);
				*(*(cmtx+1)+9+ii) = -*(cross2+ii);
			}

			dot1 = dotproduct(axis1,axis3,3);
			dot2 = dotproduct(axis2,axis3,3);

			/*
			for (ii = 0; ii < 3; ii++)
			{
				for (jj = 0; jj < 3; jj++)
				{
					*(*(cmtxH+3+ii)+3+jj) +=   *(axis1+ii)**(axis3+jj)**(lambda+leq+0)+*(axis2+ii)**(axis3+jj)**(lambda+leq+1);
					*(*(cmtxH+3+ii)+9+jj) +=  -*(axis3+ii)**(axis1+jj)**(lambda+leq+0)-*(axis3+ii)**(axis2+jj)**(lambda+leq+1);
					*(*(cmtxH+9+ii)+3+jj) +=  -*(axis1+ii)**(axis3+jj)**(lambda+leq+0)-*(axis2+ii)**(axis3+jj)**(lambda+leq+1);
					*(*(cmtxH+9+ii)+9+jj) +=   *(axis3+ii)**(axis1+jj)**(lambda+leq+0)+*(axis3+ii)**(axis2+jj)**(lambda+leq+1);
					if(ii==jj)
					{
						*(*(cmtxH+3+ii)+3+jj) -= dot1**(lambda+leq+0)+dot2**(lambda+leq+1);
						*(*(cmtxH+3+ii)+9+jj) += dot1**(lambda+leq+0)+dot2**(lambda+leq+1);
						*(*(cmtxH+9+ii)+3+jj) += dot1**(lambda+leq+0)+dot2**(lambda+leq+1);
						*(*(cmtxH+9+ii)+9+jj) -= dot1**(lambda+leq+0)+dot2**(lambda+leq+1);
					}
				}
			}
			symmetricmtx(cmtxH, 6*nnod);
			*/



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
			for (ii = 0; ii < 6 * nnod; ii++)
			{
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					if(*(*(cmtxH+ii)+jj)!=0.0)
					{
						gread(gmtx , *(constraintmain + *(loffset + ii)) + 1, *(constraintmain + *(loffset + jj)) + 1, &gdata);
						gwrite(gmtx , *(constraintmain + *(loffset + ii)) + 1, *(constraintmain + *(loffset + jj)) + 1, gdata + *(*(cmtxH+ii)+jj));
					}
				}
			}

			free(loffset);
			free(gforminit);
			free(gform);
			free(gdisp);

			free(rvct_m);
			free(rvct_s);
			freematrix(R_m,3);
			freematrix(R_s,3);

			free(initaxis1);
			free(initaxis2);
			free(initaxis3);
			free(axis1);
			free(axis2);
			free(axis3);
			free(cross1);
			free(cross2);
			freematrix(cmtx,neq);
			freematrix(cmtxH,6*nnod);
		}

		if(type == 10)/*PERIODIC BOUNDARY.*/
		{
			/*Jacobian*/
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

			gforminit = extractconstraintdisplacement(*(constraints+i-1), iform);
			gform = extractconstraintdisplacement(*(constraints+i-1), ddisp);
			gdisp = extractdeformation(gforminit, gform, nnod);



			if(nnod==3)
			{
				for (ii = 0; ii < 3; ii++)
				{
					for (jj = 0; jj < 3; jj++)
					{
						if(ii==jj)
						{
							*(*(cmtx+ii)   +jj) =  1.0;
							*(*(cmtx+ii)+6 +jj) = -1.0;
							*(*(cmtx+ii)+12+jj) =  1.0;
						}
					}
				}
			}
			if(nnod==4)
			{
				for (ii = 0; ii < 3; ii++)
				{
					for (jj = 0; jj < 3; jj++)
					{
						if(ii==jj)
						{
							*(*(cmtx+ii)   +jj) =  1.0;
							*(*(cmtx+ii)+6 +jj) = -1.0;
							*(*(cmtx+ii)+12+jj) =  1.0;
							*(*(cmtx+ii)+18+jj) =  1.0;
						}
					}
				}
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
			for (ii = 0; ii < 6 * nnod; ii++)
			{
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					if(*(*(cmtxH+ii)+jj)!=0.0)
					{
						gread(gmtx , *(constraintmain + *(loffset + ii)) + 1, *(constraintmain + *(loffset + jj)) + 1, &gdata);
						gwrite(gmtx , *(constraintmain + *(loffset + ii)) + 1, *(constraintmain + *(loffset + jj)) + 1, gdata + *(*(cmtxH+ii)+jj));
					}
				}
			}


			free(loffset);
			free(gforminit);
			free(gform);
			free(gdisp);

			freematrix(cmtx,neq);
		}


		if(type == 11)/*PERIODIC BOUNDARY.*/
		{
			/*Jacobian*/
			cmtx = (double**)malloc(neq * sizeof(double*));
			for(ii = 0 ; ii < neq ; ii++)
			{
				*(cmtx+ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6*nnod; jj++)
				{
					*(*(cmtx+ii)+jj) = 0.0;
				}
			}

			/*Hessian*/
			cmtxH = (double**)malloc(6 * nnod * sizeof(double*));
			for(ii = 0 ; ii < 6 * nnod ; ii++)/*HESSIAN OF CONSTRAINT*/
			{
				*(cmtxH+ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					*(*(cmtxH+ii)+jj) = 0.0;
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

			gforminit = extractconstraintdisplacement(*(constraints+i-1), iform);
			gform = extractconstraintdisplacement(*(constraints+i-1), ddisp);
			gdisp = extractdeformation(gforminit, gform, nnod);


			/*Rotation Vector*/
			rvct_m = (double*)malloc(3 * sizeof(double));
			rvct_s = (double*)malloc(3 * sizeof(double));
			rvct_c = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(rvct_m + ii) = *(gdisp + *(loffset + 3 + ii));
				*(rvct_s + ii) = *(gdisp + *(loffset + 9 + ii));
				*(rvct_c + ii) = *(gdisp + *(loffset + 15 + ii));
			}
			R_m = rotationmtx(rvct_m);
			R_s = rotationmtx(rvct_s);
			R_c = rotationmtx(rvct_c);


			R_st = matrixtranspose(R_s,3);
			R_cm = matrixmatrix(R_c,R_m,3);
			R_cmst = matrixmatrix(R_cm,R_st,3);
			rvct_cmst = rotationvct(R_cmst);

			H =  jacobimtx(rvct_cmst);

			HR_cmst = matrixmatrix(H,R_cmst,3);
			HR_c = matrixmatrix(H,R_c,3);

			for (ii = 0; ii < 3; ii++)
			{
				for (jj = 0; jj < 3; jj++)
				{
					*(*(cmtx+3+ii)+ 3+jj) =  *(*(HR_c+ii)+jj);
					*(*(cmtx+3+ii)+ 9+jj) = -*(*(HR_cmst+ii)+jj);
					*(*(cmtx+3+ii)+15+jj) =  *(*(H+ii)+jj);
					if(ii==jj)
					{
						*(*(cmtx+ii)   +jj) =  1.0;
						*(*(cmtx+ii)+6 +jj) = -1.0;
						*(*(cmtx+ii)+12+jj) =  1.0;
					}
				}
			}

			/*
			for (ii = 0; ii < 3; ii++)
			{
				*(*(cmtxH+3+ii)+3) += *(*(rmtx+2)+ii)*(lambda+leq+4)-*(*(rmtx+1)+ii)*(lambda+leq+5);
				*(*(cmtxH+3+ii)+4) += *(*(rmtx+0)+ii)*(lambda+leq+5)-*(*(rmtx+2)+ii)*(lambda+leq+3);
				*(*(cmtxH+3+ii)+5) += *(*(rmtx+1)+ii)*(lambda+leq+3)-*(*(rmtx+0)+ii)*(lambda+leq+4);
			}
			symmetricmtx(cmtxH, 6*nnod);
			*/


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
			for (ii = 0; ii < 6 * nnod; ii++)
			{
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					if(*(*(cmtxH+ii)+jj)!=0.0)
					{
						gread(gmtx , *(constraintmain + *(loffset + ii)) + 1, *(constraintmain + *(loffset + jj)) + 1, &gdata);
						gwrite(gmtx , *(constraintmain + *(loffset + ii)) + 1, *(constraintmain + *(loffset + jj)) + 1, gdata + *(*(cmtxH+ii)+jj));
					}
				}
			}


			free(loffset);
			free(gforminit);
			free(gform);
			free(gdisp);

			free(rvct_m);
			free(rvct_s);
			free(rvct_c);
			free(rvct_cmst);
			freematrix(R_m,3);
			freematrix(R_s,3);
			freematrix(R_c,3);
			freematrix(R_st,3);
			freematrix(R_cm,3);
			freematrix(R_cmst,3);

			freematrix(H,3);
			freematrix(HR_cmst,3);
			freematrix(HR_c,3);

			freematrix(cmtx,neq);
			freematrix(cmtxH,6*nnod);
		}


	}
	return;
}




void constraintstress(struct oconstraint* constraints, int nconstraint, long int* constraintmain,
					  double* iform, double* ddisp, double* lambda, double* fconstraint, double* constraintvct)
{
	char str[500];
	int i,j,ii,jj;
	int* loffset;
	int nnod,type,neq,leq;

	double* cvct;
	double** cmtx;
	double** cmtxH;

	double gdata;

	double* initaxis1,* initaxis2,* initaxis3;
	double* axis1,* axis2,* axis3;
	double *cross1,*cross2;
	double dot1 = 0.0;
	double dot2 = 0.0;

	double *gforminit, *gform, *gdisp;
	double* gconstraint;

	double *rvct_m, *rvct_s, *rvct_c, *rvct_cmst;
	double **R_m, **R_s, **R_c;
	double **R_st, **R_cm, **R_cmst;
	double **H, **HR_cmst, **HR_c;


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

			gforminit = extractconstraintdisplacement(*(constraints+i-1), iform);
			gform = extractconstraintdisplacement(*(constraints+i-1), ddisp);
			gdisp = extractdeformation(gforminit, gform, nnod);


			rvct_m = (double*)malloc(3 * sizeof(double));
			rvct_s = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(rvct_m + ii) = *(gdisp + *(loffset + 3 + ii));
				*(rvct_s + ii) = *(gdisp + *(loffset + 9 + ii));
			}
			R_m = rotationmtx(rvct_m);
			R_s = rotationmtx(rvct_s);

			initaxis1 = (double*)malloc(3 * sizeof(double));
			initaxis2 = (double*)malloc(3 * sizeof(double));
			initaxis3 = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(initaxis1 + ii) = (constraints+i-1)->axis[0][ii];
				*(initaxis2 + ii) = (constraints+i-1)->axis[1][ii];
				*(initaxis3 + ii) = (constraints+i-1)->axis[2][ii];
			}
			axis1 = matrixvector(R_m,initaxis1,3);
			axis2 = matrixvector(R_m,initaxis2,3);
			axis3 = matrixvector(R_s,initaxis3,3);

			dot1 = dotproduct(axis1,axis3,3);
			dot2 = dotproduct(axis2,axis3,3);
			*(cvct+0) = dot1;
			*(cvct+1) = dot2;

			cross1 = crossproduct(axis1,axis3);
			cross2 = crossproduct(axis2,axis3);
			for (ii = 0; ii < 3; ii++)
			{
				*(*(cmtx+0)+3+ii) =  *(cross1+ii);
				*(*(cmtx+0)+9+ii) = -*(cross1+ii);
				*(*(cmtx+1)+3+ii) =  *(cross2+ii);
				*(*(cmtx+1)+9+ii) = -*(cross2+ii);
			}


			gconstraint = (double*)malloc(6*nnod * sizeof(double));
			for (ii = 0; ii < 6*nnod; ii++)
			{
				*(gconstraint + ii) = 0.0;
				for (jj = 0; jj < neq; jj++)
				{
					*(gconstraint + ii) += *(lambda+leq+jj)**(*(cmtx+jj)+ii);
				}
				*(fconstraint + *(constraintmain + *(loffset + ii))) += *(gconstraint + ii);
			}

			for (ii = 0; ii < neq; ii++)
			{
				*(constraintvct+leq+ii) += *(cvct+ii);
			}


			free(loffset);
			free(gforminit);
			free(gform);
			free(gdisp);
			free(gconstraint);

			free(rvct_m);
			free(rvct_s);
			freematrix(R_m,3);
			freematrix(R_s,3);

			free(initaxis1);
			free(initaxis2);
			free(initaxis3);
			free(axis1);
			free(axis2);
			free(axis3);
			free(cross1);
			free(cross2);

			free(cvct);
			freematrix(cmtx,neq);

		}



		if(type == 10)/*PERIODIC BOUNDARY 2D.*/
		{
			/*Jacobian*/
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

			gforminit = extractconstraintdisplacement(*(constraints+i-1), iform);
			gform = extractconstraintdisplacement(*(constraints+i-1), ddisp);
			gdisp = extractdeformation(gforminit, gform, nnod);




			for (ii = 0; ii < 3; ii++)
			{
			  *(cvct+  ii) = *(gdisp+12+ii) + *(gdisp+ii) - *(gdisp+6+ii);
			}


			for (ii = 0; ii < 3; ii++)
			{
				for (jj = 0; jj < 3; jj++)
				{
					if(ii==jj)
					{
						*(*(cmtx+ii)   +jj) =  1.0;
						*(*(cmtx+ii)+6 +jj) = -1.0;
						*(*(cmtx+ii)+12+jj) =  1.0;
					}
				}
			}



			gconstraint = (double*)malloc(6*nnod * sizeof(double));
			for (ii = 0; ii < 6*nnod; ii++)
			{
				*(gconstraint + ii) = 0.0;
				for (jj = 0; jj < neq; jj++)
				{
					*(gconstraint + ii) += *(lambda+leq+jj)**(*(cmtx+jj)+ii);
				}
				*(fconstraint + *(constraintmain + *(loffset + ii))) += *(gconstraint + ii);
			}

			for (ii = 0; ii < neq; ii++)
			{
				*(constraintvct+leq+ii) += *(cvct+ii);
			}


			free(loffset);
			free(gforminit);
			free(gform);
			free(gdisp);
			free(gconstraint);

			free(rvct_m);
			free(rvct_s);
			free(rvct_c);
			free(rvct_cmst);
			freematrix(R_m,3);
			freematrix(R_s,3);
			freematrix(R_c,3);
			freematrix(R_st,3);
			freematrix(R_cm,3);
			freematrix(R_cmst,3);

			freematrix(H,3);
			freematrix(HR_cmst,3);
			freematrix(HR_c,3);

			free(cvct);
			freematrix(cmtx,neq);
		}



		if(type == 11)/*PERIODIC BOUNDARY 2D.*/
		{
			/*Jacobian*/
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

			gforminit = extractconstraintdisplacement(*(constraints+i-1), iform);
			gform = extractconstraintdisplacement(*(constraints+i-1), ddisp);
			gdisp = extractdeformation(gforminit, gform, nnod);


			/*Rotation Vector*/
			rvct_m = (double*)malloc(3 * sizeof(double));
			rvct_s = (double*)malloc(3 * sizeof(double));
			rvct_c = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(rvct_m + ii) = *(gdisp + *(loffset + 3 + ii));
				*(rvct_s + ii) = *(gdisp + *(loffset + 9 + ii));
				*(rvct_c + ii) = *(gdisp + *(loffset + 15 + ii));
			}
			R_m = rotationmtx(rvct_m);
			R_s = rotationmtx(rvct_s);
			R_c = rotationmtx(rvct_c);


			R_st = matrixtranspose(R_s,3);
			R_cm = matrixmatrix(R_c,R_m,3);
			R_cmst = matrixmatrix(R_cm,R_st,3);
			rvct_cmst = rotationvct(R_cmst);

			for (ii = 0; ii < 3; ii++)
			{
			  *(cvct+  ii) = *(gdisp+12+ii) + *(gdisp+ii) - *(gdisp+6+ii);
			  *(cvct+3+ii) = *(rvct_cmst+ii);
			}

			H =  jacobimtx(rvct_cmst);
			HR_cmst = matrixmatrix(H,R_cmst,3);
			HR_c = matrixmatrix(H,R_c,3);

			for (ii = 0; ii < 3; ii++)
			{
				for (jj = 0; jj < 3; jj++)
				{
					*(*(cmtx+3+ii)+ 3+jj) =  *(*(HR_c+ii)+jj);
					*(*(cmtx+3+ii)+ 9+jj) = -*(*(HR_cmst+ii)+jj);
					*(*(cmtx+3+ii)+15+jj) =  *(*(H+ii)+jj);
					if(ii==jj)
					{
						*(*(cmtx+ii)   +jj) =  1.0;
						*(*(cmtx+ii)+6 +jj) = -1.0;
						*(*(cmtx+ii)+12+jj) =  1.0;
					}
				}
			}



			gconstraint = (double*)malloc(6*nnod * sizeof(double));
			for (ii = 0; ii < 6*nnod; ii++)
			{
				*(gconstraint + ii) = 0.0;
				for (jj = 0; jj < neq; jj++)
				{
					*(gconstraint + ii) += *(lambda+leq+jj)**(*(cmtx+jj)+ii);
				}
				*(fconstraint + *(constraintmain + *(loffset + ii))) += *(gconstraint + ii);
			}

			for (ii = 0; ii < neq; ii++)
			{
				*(constraintvct+leq+ii) += *(cvct+ii);
			}


			free(loffset);
			free(gforminit);
			free(gform);
			free(gdisp);
			free(gconstraint);

			free(rvct_m);
			free(rvct_s);
			free(rvct_c);
			free(rvct_cmst);
			freematrix(R_m,3);
			freematrix(R_s,3);
			freematrix(R_c,3);
			freematrix(R_st,3);
			freematrix(R_cm,3);
			freematrix(R_cmst,3);

			freematrix(H,3);
			freematrix(HR_cmst,3);
			freematrix(HR_c,3);

			free(cvct);
			freematrix(cmtx,neq);
		}


	}
	return;
}





void assemconstraint_DYNA(struct oconstraint* constraints, int nconstraint, long int* constraintmain,
						  struct gcomponent* gmtx,
						  double* iform, double* lastddisp, double* ddisp, double*lastlambda, double* lambda,
						  double alphaf, int msize)
{
	char str[500];
	int i,j,ii,jj;
	int* loffset;

	int nnod,type,neq,leq;
	double gdata;

	double dot111,dot112,dot122;
	double dot211,dot212,dot222;

	double* initaxis1,* initaxis2,* initaxis3;

	double* rvct1,* rvct2;

	double** lastrmtx1,** lastrmtx2;
	double** rmtx1,** rmtx2;
	double** midrmtx1,** midrmtx2;

	double* axis1,* axis2,* axis3;
	double* midaxis1,* midaxis2,* midaxis3;

	double *cross11,*cross12;
	double *cross21,*cross22;


	double* cvct;
	double** cmtx;
	double** H;

	for (i = 1; i <= nconstraint; i++)
	{
		neq = (constraints+i-1)->neq;
		leq = (constraints+i-1)->leq;
		type = (constraints+i-1)->type;
		nnod = (constraints+i-1)->nnod;

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
			H = (double**)malloc(6 * nnod * sizeof(double*));
			for(ii = 0 ; ii < 6 * nnod ; ii++)/*HESSIAN OF CONSTRAINT*/
			{
				*(H+ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					*(*(H+ii)+jj) = 0.0;
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
				*(rvct1 + ii) = *(lastddisp + *(loffset + 3 + ii));
				*(rvct2 + ii) = *(lastddisp + *(loffset + 9 + ii));
			}
			lastrmtx1 = rotationmtx(rvct1);
			lastrmtx2 = rotationmtx(rvct2);

			for (ii = 0; ii < 3; ii++)
			{
				*(rvct1 + ii) = *(ddisp + *(loffset + 3 + ii));
				*(rvct2 + ii) = *(ddisp + *(loffset + 9 + ii));
			}
			rmtx1 = rotationmtx(rvct1);
			rmtx2 = rotationmtx(rvct2);

			midrmtx1 = midpointmtx(rmtx1, lastrmtx1, alphaf, 3);
			midrmtx1 = midpointmtx(rmtx2, lastrmtx2, alphaf, 3);


			initaxis1 = (double*)malloc(3 * sizeof(double));
			initaxis2 = (double*)malloc(3 * sizeof(double));
			initaxis3 = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(initaxis1 + ii) = (constraints+i-1)->axis[0][ii];
				*(initaxis2 + ii) = (constraints+i-1)->axis[1][ii];
				*(initaxis3 + ii) = (constraints+i-1)->axis[2][ii];
			}


			//dbgvct(initaxis1,3,3,NULL);
			//dbgvct(initaxis2,3,3,NULL);
			//dbgvct(initaxis3,3,3,NULL);


			axis1 = matrixvector(rmtx1,initaxis1,3);
			axis2 = matrixvector(rmtx1,initaxis2,3);
			axis3 = matrixvector(rmtx2,initaxis3,3);

			midaxis1 = matrixvector(midrmtx1,initaxis1,3);
			midaxis2 = matrixvector(midrmtx1,initaxis2,3);
			midaxis3 = matrixvector(midrmtx2,initaxis3,3);




			dot111 = dotproduct(axis1,midaxis3,3);
			dot112 = dotproduct(axis1,axis3,3);
			dot122 = dotproduct(midaxis1,axis3,3);

			dot211 = dotproduct(axis2,midaxis3,3);
			dot212 = dotproduct(axis2,axis3,3);
			dot222 = dotproduct(midaxis2,axis3,3);


			cross11 = crossproduct(axis1,midaxis3);
			cross12 = crossproduct(axis3,midaxis1);

			cross21 = crossproduct(axis2,midaxis3);
			cross22 = crossproduct(axis3,midaxis2);

			for (ii = 0; ii < 3; ii++)
			{
				*(*(cmtx+0)+3+ii) = *(cross11+ii);
				*(*(cmtx+0)+9+ii) = *(cross12+ii);
				*(*(cmtx+1)+3+ii) = *(cross21+ii);
				*(*(cmtx+1)+9+ii) = *(cross22+ii);
			}

			for (ii = 0; ii < 3; ii++)
			{
				for (jj = 0; jj < 3; jj++)
				{
					*(*(H+3+ii)+3+jj) +=   (1-alphaf)*           (*(axis1+ii)**(midaxis3+jj)**(lambda+leq+0)+*(axis2+ii)**(midaxis3+jj)**(lambda+leq+1));
					*(*(H+3+ii)+9+jj) +=  -(1-alphaf)*(1-alphaf)*(*(axis3+ii)*   *(axis1+jj)**(lambda+leq+0)-*(axis3+ii)*   *(axis2+jj)**(lambda+leq+1));
					*(*(H+9+ii)+3+jj) +=  -(1-alphaf)*(1-alphaf)*(*(axis1+ii)*   *(axis3+jj)**(lambda+leq+0)-*(axis2+ii)*   *(axis3+jj)**(lambda+leq+1));
					*(*(H+9+ii)+9+jj) +=   (1-alphaf)*           (*(axis3+ii)**(midaxis1+jj)**(lambda+leq+0)+*(axis3+ii)**(midaxis2+jj)**(lambda+leq+1));
					if(ii==jj)
					{
						*(*(H+3+ii)+3+jj) += -(1-alphaf)           *(dot111**(lambda+leq+0)+dot211**(lambda+leq+1));
						*(*(H+3+ii)+9+jj) +=  (1-alphaf)*(1-alphaf)*(dot112**(lambda+leq+0)+dot212**(lambda+leq+1));
						*(*(H+9+ii)+3+jj) +=  (1-alphaf)*(1-alphaf)*(dot112**(lambda+leq+0)+dot212**(lambda+leq+1));
						*(*(H+9+ii)+9+jj) += -(1-alphaf)           *(dot122**(lambda+leq+0)+dot222**(lambda+leq+1));
					}
				}
			}
			symmetricmtx(H, 6);

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
			for (ii = 0; ii < 6 * nnod; ii++)
			{
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					if(*(*(H+ii)+jj)!=0.0)
					{
						gread(gmtx , *(constraintmain + *(loffset + ii)) + 1, *(constraintmain + *(loffset + jj)) + 1, &gdata);
						gwrite(gmtx , *(constraintmain + *(loffset + ii)) + 1, *(constraintmain + *(loffset + jj)) + 1, gdata + *(*(H+ii)+jj));
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
			free(cross11);
			free(cross12);
			free(cross21);
			free(cross22);
			freematrix(cmtx,neq);
			freematrix(H,6*nnod);
		}

	}
	return;
}




void constraintstress_DYNA(struct oconstraint* constraints, int nconstraint, long int* constraintmain,
						   double* iform, double* lastddisp, double* ddisp, double* lastlambda, double* lambda,
						   double* fconstraint, double* constraintvct,
						   double alphaf)
{
	char str[500];
	int i,j,ii,jj;
	int* loffset;

	int nnod,type,neq,leq;

	double dot1 = 0.0;
	double dot2 = 0.0;


	double* initaxis1,* initaxis2,* initaxis3;

	double* rvct1,* rvct2;

	double** lastrmtx1,** lastrmtx2;
	double** rmtx1,** rmtx2;
	double** midrmtx1,** midrmtx2;

	double* axis1,* axis2,* axis3;
	double* midaxis1,* midaxis2,* midaxis3;

	double *cross11,*cross12;
	double *cross21,*cross22;


	double* gconstraint;

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
				*(rvct1 + ii) = *(lastddisp + *(loffset + 3 + ii));
				*(rvct2 + ii) = *(lastddisp + *(loffset + 9 + ii));
			}
			lastrmtx1 = rotationmtx(rvct1);
			lastrmtx2 = rotationmtx(rvct2);

			for (ii = 0; ii < 3; ii++)
			{
				*(rvct1 + ii) = *(ddisp + *(loffset + 3 + ii));
				*(rvct2 + ii) = *(ddisp + *(loffset + 9 + ii));
			}
			rmtx1 = rotationmtx(rvct1);
			rmtx2 = rotationmtx(rvct2);

			midrmtx1 = midpointmtx(rmtx1, lastrmtx1, alphaf, 3);
			midrmtx1 = midpointmtx(rmtx2, lastrmtx2, alphaf, 3);



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

			midaxis1 = matrixvector(midrmtx1,initaxis1,3);
			midaxis2 = matrixvector(midrmtx1,initaxis2,3);
			midaxis3 = matrixvector(midrmtx2,initaxis3,3);


			dot1 = dotproduct(midaxis1,midaxis3,3);
			dot2 = dotproduct(midaxis2,midaxis3,3);
			*(cvct+0) = dot1;
			*(cvct+1) = dot2;

			cross11 = crossproduct(axis1,midaxis3);
			cross12 = crossproduct(axis3,midaxis1);

			cross21 = crossproduct(axis2,midaxis3);
			cross22 = crossproduct(axis3,midaxis2);

			for (ii = 0; ii < 3; ii++)
			{
				*(*(cmtx+0)+3+ii) = *(cross11+ii);
				*(*(cmtx+0)+9+ii) = *(cross12+ii);
				*(*(cmtx+1)+3+ii) = *(cross21+ii);
				*(*(cmtx+1)+9+ii) = *(cross22+ii);
			}


			gconstraint = (double*)malloc(6*nnod * sizeof(double));
			for (ii = 0; ii < 6*nnod; ii++)
			{
				*(gconstraint + ii) = 0.0;
				for (jj = 0; jj < neq; jj++)
				{
					*(gconstraint + ii) += ((1-alphaf)**(lambda+leq+jj)+alphaf**(lastlambda+leq+jj))**(*(cmtx+jj)+ii);
				}
				*(fconstraint + *(constraintmain + *(loffset + ii))) += *(gconstraint + ii);
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
			free(cross11);
			free(cross12);
			free(cross21);
			free(cross22);
			freematrix(cmtx,neq);
			free(cvct);
			free(gconstraint);
		}

	}
	return;
}


