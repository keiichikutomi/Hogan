
void assemelem(struct owire* elems, struct memoryelem* melem, int nelem, long int* constraintmain,
			   struct gcomponent* mmtx, struct gcomponent* gmtx,
			   double* iform, double* ddisp, double* finternal, double* fexternal)
{
	struct owire elem;
	int i,j,ii,jj;
	int nnod;
	double* gforminit, * eforminit;                      /*DEFORMED COORDINATION OF ELEMENT*/
	double* gform, * eform;                      /*INITIAL COORDINATION OF ELEMENT*/
	double* edisp;
	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	//double* gexternal, * eexternal;                          /*EXTERNAL FORCE OF ELEMENT*/
	double** Me,** Ke,** Kt,** DBe,** B,** drccos,** drccosinit;                           /*MATRIX*/
	double** T,** Tt, **HPT,**TtPtHt;
	int* loffset;
	double area;


	for (i = 1; i <= nelem; i++)
	{
		inputelem(elems,melem,i-1,&elem);
		nnod=elem.nnod;
		elem.sect=(elems+i-1)->sect;
		loffset = (int*)malloc(6 * nnod * sizeof(int));
		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				*(loffset + (6 * ii + jj)) = 6 * (elem.node[ii]->loff) + jj;
			}
		}

		/*INITIAL CONFIGURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(iform, elem.node[ii]);
		}
		drccosinit = directioncosine(elem.node[0]->d[0],
							elem.node[0]->d[1],
							elem.node[0]->d[2],
							elem.node[1]->d[0],
							elem.node[1]->d[1],
							elem.node[1]->d[2],
							elem.cangle);
		gforminit = extractdisplacement(elem, iform);
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);

		Ke = assememtx(elem);
		Ke = modifyhinge(elem,Ke);             /*MODIFY MATRIX.*/

		/*DEFORMED CONFIDURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, elem.node[ii]);
		}
		gform = extractdisplacement(elem, ddisp);
		drccos = updatedrccos(drccosinit, gforminit, gform);
		eform = extractlocalcoord(gform,drccos,nnod);

		T = transmatrixIII(drccos, nnod);									/*[T].*/
		Tt = matrixtranspose(T, 6 * nnod);                    				/*[Tt].*/

		edisp = extractdeformation(eforminit, eform, nnod);                	/*{Ue}*/
		einternal = matrixvector(Ke, edisp, 6 * nnod);          			/*{Fe}=[Ke]{Ue}.*/

		ginternal = (double*)malloc(6 * nnod * sizeof(double));
		HPT = (double**)malloc(6 * nnod * sizeof(double*));
		for (ii = 0; ii <6 * nnod; ii++)
		{
			*(HPT + ii) = (double*)malloc(6 * nnod * sizeof(double));
		}
		Kt = assemtmtxCR(Ke, eform, edisp, einternal, ginternal, T, HPT, nnod);	/*TANGENTIAL MATRIX[Kt].*/
		symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/

		if(gmtx!=NULL)assemgstiffnesswithDOFelimination(gmtx, Kt, &elem, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/



		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				if(finternal!=NULL)
				{
					*(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				}
				(elems+i-1)->stress[ii][jj]=*(einternal+6*ii+jj);
				if(melem!=NULL)
				{
					(melem+i-1)->stress[ii][jj]=*(einternal+6*ii+jj);
				}
			}
		}




		freematrix(drccosinit, 3);
		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(Ke, 6 * nnod);
		freematrix(Kt, 6 * nnod);

		free(einternal);
		free(ginternal);

		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);
		free(loffset);
	}
	return;
}

void assemshell(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain,
				struct gcomponent* mmtx, struct gcomponent* gmtx,
				double* iform, double* ddisp, double* finternal, double* fexternal)
{
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	int ngp;
	int nstress;

	double* gforminit, * eforminit;                      /*DEFORMED COORDINATION OF ELEMENT*/
	double* gform, * eform;                      /*INITIAL COORDINATION OF ELEMENT*/
	double* edisp;

	double** drccos,** drccosinit;                           /*MATRIX*/
	double** T,** Tt;

	double** estress, ** estrain;
	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	double* gexternal, * eexternal;                          /*EXTERNAL FORCE OF ELEMENT*/

	double*** C, ***B;
	double** M,** Kp, **Ke, ** Kt;
	double** HPT,** TtPtHt;

	int* loffset;
	double area;


	for (i = 1; i <= nshell; i++)
	{

		inputshell(shells, mshell, i - 1, &shell);

		nnod = shell.nnod;
		ngp = shell.ngp;
		nstress = shell.nstress;

		loffset = (int*)malloc(6 * nnod * sizeof(int));
		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				*(loffset + (6 * ii + jj)) = 6 * (shell.node[ii]->loff) + jj;
			}
		}

		/*INITIAL CONFIGURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(iform, shell.node[ii]);
		}
		drccosinit = shelldrccos(shell);
		gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);        	/*{Xe}*/

		C = elasticCshell(shell);
		B = assemshellshape(shell, drccosinit);
		//Ke = assemshellemtx(shell);                                      /*[Ke]*/
		M = assemshellmmtx(shell, drccosinit);         				       /*[Me]*/

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
		Tt = matrixtranspose(T, 6 * nnod);                  				/*[Tt].*/

		HPT = transmatrixHPT(eform, edisp, T, nnod);
		TtPtHt = matrixtranspose(HPT, 6 * nnod);

		assemshellestrain(&shell, B, edisp);
		assemshellestress(&shell, C);
		einternal = assemshelleinternal(&shell, B);

		if(finternal!=NULL)
		{
			//einternal = matrixvector(Ke, edisp, 6*nnod);

			ginternal = matrixvector(TtPtHt, einternal, 6 * nnod);
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				}
			}

			free(ginternal);
		}
		if(fexternal!=NULL)
		{
			eexternal = assemshellpvct(shell, drccos);
			gexternal = matrixvector(Tt, eexternal, 6 * nnod);
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(fexternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(gexternal + 6 * ii + jj);
				}
			}
			free(eexternal);
			free(gexternal);
		}
		if(gmtx!=NULL)
		{

			Kp = assemshellpmtx(shell,C,B);
			Kp = transformationIII(Kp, HPT, 6*nnod);/*[Ke]=[Tt][Pt][Ht][K][H][P][T]*/
			Kt = assemgmtxCR(eform, edisp, einternal, NULL, T, NULL, nnod);/*[Kg]=[Kgr]+[Kgp]+[Kgm]*/
			for (ii = 0; ii < 6*nnod; ii++)
			{
				for (jj = 0; jj < 6*nnod; jj++)
				{
					*(*(Kt + ii) + jj) += *(*(Kp + ii) + jj);/*[Kt]=[Ke]+[Kg]*/
				}
			}
			symmetricmtx(Kt, 6 * nnod);											   /*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
			assemgstiffnessIIwithDOFelimination(gmtx, Kt, &shell, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/
			freematrix(Kp, 6 * nnod);
			freematrix(Kt, 6 * nnod);
		}
		if(mmtx!=NULL)
		{
			M = transformationIII(M, T, 6*nnod);
			assemgstiffnessIIwithDOFelimination(mmtx, M, &shell, constraintmain);
		}

		//volume += shellvolume(shell, drccos, area);                   		/*VOLUME*/
		outputshell(shells, i - 1, &shell);

		for(ii=0;ii<ngp;ii++)freematrix(*(C+ii), nstress);
		free(C);
		for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
		free(B);

        freematrix(M, 6 * nnod);

		free(loffset);
		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);
		freematrix(drccos, 3);
		freematrix(drccosinit, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(TtPtHt, 6 * nnod);

		free(einternal);
	}

	return;
}


void assemshellEx(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain, struct oconf* confs,
				  std::vector<Triplet>& Mtriplet, std::vector<Triplet>& Ktriplet,
				  double* iform, double* ddisp, double* finternal, double* fexternal)
{
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	int ngp;
	int nstress;

	double* gforminit, * eforminit;                      /*DEFORMED COORDINATION OF ELEMENT*/
	double* gform, * eform;                      /*INITIAL COORDINATION OF ELEMENT*/
	double* edisp;

	double** drccos,** drccosinit;                           /*MATRIX*/
	double** T,** Tt;

	double** estress, ** estrain;
	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	double* gexternal, * eexternal;                          /*EXTERNAL FORCE OF ELEMENT*/

	double*** C, ***B;
	double** M,** K, **Ke, ** Kt;
	double** HPT,** TtPtHt;

	int* loffset;
	double area;


	for (i = 1; i <= nshell; i++)
	{

		inputshell(shells, mshell, i - 1, &shell);

		nnod = shell.nnod;
		ngp = shell.ngp;
		nstress = shell.nstress;

		loffset = (int*)malloc(6 * nnod * sizeof(int));
		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				*(loffset + (6 * ii + jj)) = 6 * (shell.node[ii]->loff) + jj;
			}
		}

		/*INITIAL CONFIGURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(iform, shell.node[ii]);
		}
		drccosinit = shelldrccos(shell);
		gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);        	/*{Xe}*/

		C = elasticCshell(shell);
		B = assemshellshape(shell, drccosinit);
		//K = assemshellemtx(shell);

		M = assemshellmmtx(shell, drccosinit);         				/*[Me]*/

		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell);
		gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
		eform = extractlocalcoord(gform,drccos,nnod); 			       	    /*{Xe+Ue}*/

		T = transmatrixIII(drccos, nnod);         							/*[T].*/
		Tt = matrixtranspose(T, 6 * nnod);                  				/*[Tt].*/

		edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/

		assemshellestrain(&shell, B, edisp);
		assemshellestress(&shell, C);
		einternal = assemshelleinternal(&shell, B);

		K = assemshellpmtx(shell,C,B);

		//einternal = matrixvector(K, edisp, 6*nnod);
		eexternal = assemshellpvct(shell, drccos);                			/*{Pe}.*/
		//volume += shellvolume(shell, drccos, area);                   	/*VOLUME*/

		ginternal = (double*)malloc(6 * nnod * sizeof(double));
		HPT = (double**)malloc(6 * nnod * sizeof(double*));
		for (ii = 0; ii < 6 * nnod; ii++)
		{
			*(HPT + ii) = (double*)malloc(6 * nnod * sizeof(double));
		}
		Kt = assemtmtxCR(K, eform, edisp, einternal, ginternal, T, HPT, nnod);	/*TANGENTIAL MATRIX[Kt].*/
		symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/

		M = transformationIII(M, T, 6*nnod);

		assemtripletshell(Ktriplet, Kt, &shell, constraintmain, confs);
		assemtripletshell(Mtriplet, M, &shell, constraintmain, confs);


		gexternal = matrixvector(Tt, eexternal, 6 * nnod);  /*GLOBAL EXTERNAL FORCE{Pg}.*/

		outputshell(shells, i - 1, &shell);

		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				if(finternal!=NULL) *(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				if(fexternal!=NULL) *(fexternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(gexternal + 6 * ii + jj);
			}
		}

		freematrix(drccos, 3);
		freematrix(drccosinit, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);

		freematrix(HPT, 6 * nnod);
		freematrix(K, 6 * nnod);
		freematrix(Kt, 6 * nnod);

		for(ii=0;ii<ngp;ii++)freematrix(*(C+ii), nstress);
		free(C);
		for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
		free(B);

		freematrix(M, 6 * nnod);

		free(einternal);
		free(ginternal);

		free(eexternal);
		free(gexternal);

		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);

		free(loffset);
	}
	return;
}



