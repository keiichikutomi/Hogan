struct oconstraint{long int code,loff;
				   int type;
				   int nequation;
				   struct onode *(node[2]);
				   double taxis1[3];
				   double taxis2[3];
				   double raxis1[3];
				   double raxis2[3];
				  };


/*REVOLUTE JOINT    : 1 DIRECTION ROTATION, NO TRANSITION.          1 VECTOR INPUT.*/
/*SPHERICAL JOINT   : FREE ROTATION, NO TRANSITION.      		    0 VECTOR INPUT.*/
/*PRISMATIC JOINT   : 1 DIRECTION TRANSITION, NO ROTATION.		    1 VECTOR INPUT.*/
/*CYLINDRICAL JOINT : 1 DIRECTION ROTATION, 1 DIRECTION TRANSITION. 2 VECTOR INPUT.*/
/*UNIVERSAL JOINT   : 2 DIRECTION ROTATION, NO TRANSITION.          2 VECTOR INPUT.*/



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

	double** drccosinit;
	double** drccos,** T,** HPT;

	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	double* gexternal, * eexternal;                        /*EXTERNAL FORCE OF ELEMENT*/




	for (i = 1; i <= nconstraint; i++)
	{
		if(constraint->type == 1)
		{
			loffset = (int*)malloc(12 * sizeof(int));
			for (ii = 0; ii < 2; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(loffset + (6 * ii + jj)) = 6 * (constraint.node[ii]->loff) + jj;
				}
			}

			/*INITIAL CONFIGURATION*/
			for (ii = 0; ii < 2; ii++)
			{
				inputnode(iform, (constraints+i-1)->node[ii]);
			}


		}




		drccosinit = shelldrccos(shell);
		gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);        	/*{Xe}*/



		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell);
		gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
		eform = extractlocalcoord(gform,drccos,nnod); 			       	    /*{Xe+Ue}*/

		edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/

		T = transmatrixIII(drccos, nnod);         							/*[T].*/
		HPT = transmatrixHPT(eform, edisp, T, nnod);



		assemgstiffnessIIwithDOFelimination(gmtx, Kt, &shell, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/





		free(loffset);
		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);
		free(einternal);
		freematrix(drccos, 3);
		freematrix(drccosinit, 3);
		freematrix(T, 6 * nnod);
		freematrix(HPT, 6 * nnod);
	}
	return;
}
