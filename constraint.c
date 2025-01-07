


/*TYPE |      JOINT        : DESCRIPTION                                                  */
/*  1  | REVOLUTE JOINT    : 1 DIRECTION ROTATION, NO TRANSITION.          1 VECTOR INPUT.*/
/*  2  | SPHERICAL JOINT   : FREE ROTATION, NO TRANSITION.      		   0 VECTOR INPUT.*/
/*  3  | PRISMATIC JOINT   : 1 DIRECTION TRANSITION, NO ROTATION.		   1 VECTOR INPUT.*/
/*  4  | CYLINDRICAL JOINT : 1 DIRECTION ROTATION, 1 DIRECTION TRANSITION. 2 VECTOR INPUT.*/
/*  5  | UNIVERSAL JOINT   : 2 DIRECTION ROTATION, NO TRANSITION.          2 VECTOR INPUT.*/


void assemconstraint(struct oconstraint* constraints, int nconstraint, long int* constraintmain,
					 struct gcomponent* mmtx, struct gcomponent* gmtx,
					 double* iform, double* ddisp)
{
	char str[500];
	int i,j,ii,jj;
	int nnod;
	int* loffset;

	double* gforminit, * eforminit;
	double* gform, * eform, * edisp;

	double* ginternal, * einternal;
	double* gexternal, * eexternal;




	for (i = 1; i <= nconstraint; i++)
	{
		if(constraint->type == 1)
		{
			loffset = (int*)malloc(12 * sizeof(int));
			for (ii = 0; ii < 2; ii++)
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

			rvct1 = (double*)malloc(3 * sizeof(double));
			rvct2 = (double*)malloc(3 * sizeof(double));
			for (ii = 0; ii < 3; ii++)
			{
				*(rvct1 + ii) = raxis[];
				*(rvct2 + ii) = *(ddisp + *(loffset + 9 + ii));
			}

			axis1 = matrixvector(rmtx1,3);
			axis2 = matrixvector(rmtx2,3);





		}




		assemgstiffnessIIwithDOFelimination(gmtx, Kt, &shell, constraintmain);

		free(loffset);
		free(rvct1);
		free(rvct2);
		freematrix(rmtx1, 3);
		freematrix(rmtx2, 3);
	}
	return;
}

