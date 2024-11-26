#if 0
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
	char str[500];
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	int ngp;
	int nstress;
	int* loffset;

	double* gforminit, * eforminit;
	double* gform, * eform, * edisp;

	double** drccosinit;
	double** drccos, ** T, ** Tt, ** HPT,** TtPtHt;;

	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	double* gexternal, * eexternal;                          /*EXTERNAL FORCE OF ELEMENT*/

	double*** C, ***B;
	double** M,** Kp, ** Kt;




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

		C = shellC(shell);
		B = shellB(shell);
		//Kp = assemshellemtx(shell);                                      /*[Ke]*/
		M = assemshellmmtx(shell);         				       /*[Me]*/

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
		//einternal = matrixvector(Kp, edisp, 6*nnod);

		if(finternal!=NULL)
		{

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

		C = shellC(shell);
		B = shellB(shell);
		//K = assemshellemtx(shell);
		M = assemshellmmtx(shell);         				/*[Me]*/

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
#endif





















/*
##########################################################################################################################################
######    ########        ####        ##            ##      ##########                  ##      ######            ##      ##########
######    ########  ####    ##  ####    ##  ######  ####    ##########    ####  ######  ####  ##########  ######  ####    ##########    ##
######    ######    ########    ##########  ############      ######      ####  ############  ##########  ############      ######      ##
####  ##    ######  ##########  ##########    ##  ######  ##  ######  ##  ####    ##  ######  ##########    ##  ######  ##  ######  ##  ##
####  ####  ########      ######      ####        ######  ##  ######  ##  ####        ######  ##########        ######  ##  ######  ##  ##
####        ############    ########    ##  ####  ######  ####  ##  ####  ####  ####  ######  ##########  ####  ######  ####  ##  ####  ##
##  ######    ############  ##########  ##  ############  ####  ##  ####  ####  ############  ##########  ############  ####  ##  ####  ##
##  ########  ####  ####    ##  ####    ##  ######  ####  ####    ######  ####  ######  ####  ######  ##  ######  ####  ####    ######  ##
	  ####      ##        ####        ##            ##      ####  ####                  ##                        ##      ####  ####
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
*/
void assemelem(struct owire* elems, struct memoryelem* melem, int nelem, long int* constraintmain,
			   struct gcomponent* mmtx, struct gcomponent* gmtx,
			   double* iform, double* ddisp)
{
	struct owire elem;
	int i,j,ii,jj;
	int nnod;
	int* loffset;

	double* gforminit, * eforminit;
	double* gform, * eform, * edisp;

	double* ginternal, * einternal;
	//double* gexternal, * eexternal;

	double** drccosinit;
	double** drccos, ** T, ** Tt, **HPT, **TtPtHt;

	double** M,** Ke,** Kp,** Kt;



	for (i = 1; i <= nelem; i++)
	{
		inputelem(elems,melem,i-1,&elem);
		nnod=elem.nnod;
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

		edisp = extractdeformation(eforminit, eform, nnod);                	/*{Ue}*/

		T = transmatrixIII(drccos, nnod);									/*[T].*/
		//Tt = matrixtranspose(T, 6 * nnod);                    				/*[Tt].*/
		HPT = transmatrixHPT(eform, edisp, T, nnod);
		//TtPtHt = matrixtranspose(HPT, 6 * nnod);

		einternal = matrixvector(Ke, edisp, 6 * nnod);          			/*{Fe}=[Ke]{Ue}.*/

		Kt = assemtmtxCR(Ke, eform, edisp, einternal, NULL, T, NULL, nnod);	/*TANGENTIAL MATRIX[Kt].*/
		symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/

		assemgstiffnesswithDOFelimination(gmtx, Kt, &elem, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/


		free(loffset);
		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);

		free(einternal);
		//free(ginternal);

		freematrix(drccosinit, 3);

		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(TtPtHt, 6 * nnod);

		//freematrix(M, 6 * nnod);
		freematrix(Ke, 6 * nnod);
		//freematrix(Kp, 6 * nnod);
		freematrix(Kt, 6 * nnod);

	}
	return;
}

/*
############################################################################################################################################
            ##      ######            ##      ##########      ##        ##              ##        ####            ####        ####        ##
##  ######  ####  ##########  ######  ####    ##########    ####  ####      ####  ####  ##  ######  ####  ######  ####  ####    ##  ####
##  ############  ##########  ############      ######      ##    ##############  ########  ######  ####  ##########    ########    ########
##    ##  ######  ##########    ##  ######  ##  ######  ##  ####  ##############  ########  ######  ####    ##  ######  ##########  ########
##        ######  ##########        ######  ##  ######  ##  ######      ########  ########        ######        ########      ######      ##
##  ####  ######  ##########  ####  ######  ####  ##  ####  ##########    ######  ########  ####  ######  ####  ############    ########
##  ############  ##########  ############  ####  ##  ####  ############  ######  ########  ####    ####  ####################  ##########
##  ######  ####  ######  ##  ######  ####  ####    ######  ####  ####    ######  ########  ######  ####  ######  ####  ####    ##  ####
            ##                        ##      ####  ####      ##        ######      ####      ####                ####        ####        ##
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
*/
void elemstress(struct owire* elems, struct memoryelem* melem, int nelem, long int* constraintmain,
				double* iform, double* ddisp, double* finternal, double* fexternal)
{
	struct owire elem;
	int i,j,ii,jj;
	int nnod;
	int* loffset;

	double* gforminit, * eforminit;
	double* gform, * eform, * edisp;

	double* ginternal, * einternal;
	//double* gexternal, * eexternal;

	double** drccosinit;
	double** drccos,** T,** Tt,**HPT,**TtPtHt;

	double** M,** Ke,** Kp,** Kt;

	for (i = 1; i <= nelem; i++)
	{
		inputelem(elems,melem,i-1,&elem);
		nnod=elem.nnod;
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

		edisp = extractdeformation(eforminit, eform, nnod);                	/*{Ue}*/

		T = transmatrixIII(drccos, nnod);									/*[T].*/
		Tt = matrixtranspose(T, 6 * nnod);                    				/*[Tt].*/
		HPT = transmatrixHPT(eform, edisp, T, nnod);
		TtPtHt = matrixtranspose(HPT, 6 * nnod);

		einternal = matrixvector(Ke, edisp, 6 * nnod);          			/*{Fe}=[Ke]{Ue}.*/
		ginternal = matrixvector(TtPtHt, einternal, 6 * nnod);

		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				*(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
			}
		}

		free(loffset);
		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);

		free(einternal);
		free(ginternal);

		freematrix(drccosinit, 3);

		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(TtPtHt, 6 * nnod);

		//freematrix(M, 6 * nnod);
		freematrix(Ke, 6 * nnod);
		//freematrix(Kp, 6 * nnod);
		//freematrix(Kt, 6 * nnod);
	}
	return;
}
/*
################################################################################################################################################
######    ########        ####        ##            ##      ##########      ##        ##      ######                  ##      ######      ######
######    ########  ####    ##  ####    ##  ######  ####    ##########    ####  ####    ##  ##########  ####  ######  ####  ##########  ########
######    ######    ########    ##########  ############      ######      ##    ##########  ##########  ####  ############  ##########  ########
####  ##    ######  ##########  ##########    ##  ######  ##  ######  ##  ####  ##########    ########  ####    ##  ######  ##########  ########
####  ####  ########      ######      ####        ######  ##  ######  ##  ######      ####              ####        ######  ##########  ########
####        ############    ########    ##  ####  ######  ####  ##  ####  ##########    ##  ##########  ####  ####  ######  ##########  ########
##  ######    ############  ##########  ##  ############  ####  ##  ####  ############  ##  ##########  ####  ############  ##########  ########
##  ########  ####  ####    ##  ####    ##  ######  ####  ####    ######  ####  ####    ##  ##########  ####  ######  ####  ######  ##  ######
      ####      ##        ####        ##            ##      ####  ####      ##        ##      ######                  ##
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
*/

void assemshell(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain,
				struct gcomponent* mmtx, struct gcomponent* gmtx,
				double* iform, double* ddisp)
{
	char str[500];
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	int ngp;
	int nstress;
	int* loffset;

	double* gforminit, * eforminit;
	double* gform, * eform, * edisp;

	double** drccosinit;
	double** drccos,** T,** Tt,** HPT,** TtPtHt;

	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	double* gexternal, * eexternal;                        /*EXTERNAL FORCE OF ELEMENT*/

	double*** C, ***B;
	double** M, ** Kp, ** Kt;


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

		//C = shellC(shell);
		B = shellB(shell);
		//Kp = assemshellemtx(shell);
		M = assemshellmmtx(shell);

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
		//Tt = matrixtranspose(T, 6 * nnod);
		HPT = transmatrixHPT(eform, edisp, T, nnod);
		//TtPtHt = matrixtranspose(HPT, 6 * nnod);

		//assemshellestrain(&shell, B, edisp);
		//assemshellestress(&shell, C);

		inputshell(shells, NULL, i - 1, &shell);
		C = shellCconsistentilyushin(shell);
		einternal = assemshelleinternal(&shell, B);

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

		if(mmtx!=NULL)
		{
		  M = transformationIII(M, T, 6*nnod);
		  assemgstiffnessIIwithDOFelimination(mmtx, M, &shell, constraintmain);
		}


		for(ii=0;ii<ngp;ii++)freematrix(*(C+ii), nstress);
		free(C);
		for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
		free(B);

		freematrix(M, 6 * nnod);
		freematrix(Kp, 6 * nnod);
		freematrix(Kt, 6 * nnod);

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
		//freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		//freematrix(TtPtHt, 6 * nnod);
	}

	return;
}

/*
##################################################################################################################################################
##        ##      ######                  ##      ######      ########        ##              ##        ####            ####        ####        ##
##  ####    ##  ##########  ####  ######  ####  ##########  ##########  ####      ####  ####  ##  ######  ####  ######  ####  ####    ##  ####
    ##########  ##########  ####  ############  ##########  ########    ##############  ########  ######  ####  ##########    ########    ########
##  ##########    ########  ####    ##  ######  ##########  ##########  ##############  ########  ######  ####    ##  ######  ##########  ########
####      ####              ####        ######  ##########  ############      ########  ########        ######        ########      ######      ##
########    ##  ##########  ####  ####  ######  ##########  ################    ######  ########  ####  ######  ####  ############    ########
##########  ##  ##########  ####  ############  ##########  ##################  ######  ########  ####    ####  ####################  ##########
##  ####    ##  ##########  ####  ######  ####  ######  ##  ######  ##  ####    ######  ########  ######  ####  ######  ####  ####    ##  ####
##        ##      ######                  ##                        ##        ######      ####      ####                ####        ####        ##
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
*/
void shellstress(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain,
				 double* iform, double* ddisp, double* finternal, double* fexternal)
{
	char str[500];
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	int ngp;
	int nstress;
	int* loffset;

	double* gforminit, * eforminit, * gform, * eform, * edisp;

	double** drccosinit;
	double** drccos,** T,** Tt,** HPT,** TtPtHt;

	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	double* gexternal, * eexternal;                          /*EXTERNAL FORCE OF ELEMENT*/

	double*** C, ***B;
	double** M,** Kp;

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

		C = shellC(shell);
		B = shellB(shell);//Kp = assemshellemtx(shell);
		M = assemshellmmtx(shell);

		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell);
		gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
		eform = extractlocalcoord(gform,drccos,nnod); 			       	    /*{Xe+Ue}*/

		edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/

		T = transmatrixIII(drccos, nnod);
		Tt = matrixtranspose(T, 6 * nnod);
		HPT = transmatrixHPT(eform, edisp, T, nnod);
		TtPtHt = matrixtranspose(HPT, 6 * nnod);

		assemshellestrain(&shell, B, edisp);
		assemshellestress(&shell, C);
		einternal = assemshelleinternal(&shell, B);//einternal = matrixvector(Kp, edisp, 6*nnod);
		eexternal = assemshellpvct(shell, drccos);

		ginternal = matrixvector(TtPtHt, einternal, 6 * nnod);
		gexternal = matrixvector(Tt, eexternal, 6 * nnod);

		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				*(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				*(fexternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(gexternal + 6 * ii + jj);
			}
		}

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
		free(einternal);
		free(ginternal);
		free(eexternal);
		free(gexternal);
		freematrix(drccosinit, 3);
		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(TtPtHt, 6 * nnod);
	}

	return;
}
















/*
##########################################################################################################################################################################################
######  ########      ####      ##          ##      ########    ##      ##      ####    ##          ##      ####      ############          ####    ####          ####      ######  ######
####  ########  ######    ######  ##  ####  ####    ######    ##  ######  ##  ########  ####  ####  ####  ########  ################  ####    ####  ######  ####    ####  ######  ########
########  ######  ########  ########  ##  ######  ########    ####  ########  ######    ####  ##  ######  ########  ################  ######  ######  ##  ######  ######  ##########  ####
####      ######      ####      ####      ######  ##  ######  ####      ####            ####      ######  ########  ################  ######  ########  ########  ##  ##  ######      ####
##  ####  ############  ########  ##  ##  ######  ##  ##  ##  ##########  ##  ########  ####  ##  ######  ########  ################  ######  ########  ########  ####    ####  ####  ####
##  ######  ##  ######    ######  ##  ##########  ####    ##  ##  ######  ##  ########  ####  ##########  ####  ##  ####  ##########  ####    ########  ########  ####    ####  ######  ##
    ####      ##      ####      ##          ##      ##  ####    ##      ##      ####    ##          ##                    ########          ########      ####      ####  ##    ####
##########################################################################################################################        ########################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
*/
void assemshell_DYNA(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain,
					 struct gcomponent* gmtx,
					 double* iform, double* lastddisp, double* ddisp, double* lapddisp,
					 double* lastudd_m, double* udd_m,
					 double* lastud_m, double* ud_m,
					 double* finertial, double* fdamping, double* finternal, double* fexternal,
					 double alpham, double alphaf, double xi, double beta, double ddt)
{
	char str[500];
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	int ngp;
	int nstress;
	int* loffset;

	double* gforminit, * eforminit;
	double* lastgform, * lasteform, * lastedisp, * lastgacc_m;
	double* gform, * eform, * edisp, * gacc_m;

	double* lasteinternal, *lasteexternal, * lastginertial_m, * lastginertial;
	double* einternal, * eexternal, * ginertial_m, * ginertial;
	double* mideinternal, * mideexternal;
	double* midginternal, * midgexternal, * midginertial;

	double** drccosinit;
	double** lastdrccos,** lastT,** lastTt,** lastHPT,** lastTtPtHt,** lastRt;
	double** drccos,** T,** Tt,** HPT,** TtPtHt,** R;
	double** midT, ** midTt, ** midHPT, ** midTtPtHt;

	double** lapH;
	double* lapgform;

	double*** C, ***B;
	double** M, ** Kp, ** Keff;

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

		//C = shellC(shell);
		B = shellB(shell);
		//Kp = assemshellemtx(shell);                                      /*[Ke]*/
		M = assemshellmmtx(shell);         				                   /*[Me]*/

		/*DEFORMED CONFIGURATION OF LAST LAP.*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(lastddisp, shell.node[ii]);
		}
		lastdrccos = shelldrccos(shell);
		lastgform = extractshelldisplacement(shell, lastddisp);             /*{Xg+Ug}*/
		lasteform = extractlocalcoord(lastgform,lastdrccos,nnod);           /*{Xe+Ue}*/

		lastedisp = extractdeformation(eforminit, lasteform, nnod);           		/*{Ue}*/

		lastT = transmatrixIII(lastdrccos, nnod);         					/*[T]*/
		lastHPT = transmatrixHPT(lasteform, lastedisp, lastT, nnod);
		lastRt = pullbackmtx(lastgform, nnod);

		inputshell(shells, mshell, i - 1, &shell);
		//C = shellCconsistentilyushin(shell);
		lasteinternal = assemshelleinternal(&shell, B);//mshell ‚©‚ç

		lastgacc_m = extractshelldisplacement(shell, lastudd_m);
		lastginertial_m = matrixvector(M, lastgacc_m, 6 * nnod);
		lastginertial = pushforward(lastgform, lastginertial_m, nnod);


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
		R = pushforwardmtx(gform,nnod);

		//assemshellestrain(&shell, B, edisp);
		//assemshellestress(&shell, C);
		inputshell(shells, NULL, i - 1, &shell);
		C = shellCconsistentilyushin(shell);
		einternal = assemshelleinternal(&shell, B);

		gacc_m = extractshelldisplacement(shell, udd_m);
		ginertial_m = matrixvector(M, gacc_m, 6 * nnod);
		ginertial = pushforward(gform, ginertial_m, nnod);

		Kp = assemshellpmtx(shell,C,B);


		lapgform = extractshelldisplacement(shell, lapddisp);                     /*{Xg+Ug}*/
		lapH = blockjacobimtx(lapgform, NULL, NULL, nnod);


			/*(24):MID-POINT TRANSFORMATION MATRIX.*/
			midHPT = midpointmtx(HPT, lastHPT, alphaf, 6 * nnod);
			midTtPtHt = matrixtranspose(midHPT, 6 * nnod);

			midT = midpointmtx(T, lastT, alphaf, 6 * nnod);
			midTt = matrixtranspose(midT, 6 * nnod);

			/*(21)&(22)&(23):MID-POINT FORCE VECTOR.*/
			mideinternal = midpointvct(einternal, lasteinternal, alphaf-xi, 6*nnod);
			midginternal = matrixvector(midTtPtHt, mideinternal, 6 * nnod);

			midginertial = midpointvct(ginertial, lastginertial, alpham, 6*nnod);

			mideexternal = midpointvct(eexternal, lasteexternal, alphaf, 6*nnod);
			midgexternal = matrixvector(midTt, mideexternal, 6 * nnod);




		double* gvel_m, * gmomentum_m;
		double masstotal,massdiag;
		double momentumlinear,momentumangular;

		Keff = assemtmtxCR_DYNA(eform, edisp, mideinternal, T, Kp, midTtPtHt, HPT,
							  ginertial, M, R, lastRt, lapH,
							  alphaf, alpham, xi, beta, ddt, nnod);

		symmetricmtx(Keff, 6*nnod);
		assemgstiffnessIIwithDOFelimination(gmtx, Keff, &shell, constraintmain);

		for(ii=0;ii<ngp;ii++)freematrix(*(C+ii), nstress);
		free(C);
		for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
		free(B);

		freematrix(M, 6 * nnod);
		freematrix(Kp, 6 * nnod);
		freematrix(Keff, 6 * nnod);

        free(loffset);

		free(eforminit);
		free(gforminit);

		free(lasteform);
		free(lastgform);
		free(lastedisp);
		free(lastgacc_m);

		free(eform);
		free(gform);
		free(edisp);
		free(gacc_m);

		free(lasteinternal);
		free(lastginertial_m);
		free(lastginertial);

		free(einternal);
		free(ginertial_m);
		free(ginertial);

		free(mideinternal);
		free(midginternal);
		free(midginertial);
		free(mideexternal);
		free(midgexternal);

		freematrix(drccosinit, 3);

		freematrix(lastdrccos, 3);
		freematrix(lastT, 6 * nnod);
		freematrix(lastTt, 6 * nnod);
		freematrix(lastHPT, 6 * nnod);
		freematrix(lastTtPtHt, 6 * nnod);
		freematrix(lastRt, 6 * nnod);

		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(TtPtHt, 6 * nnod);
		freematrix(R, 6 * nnod);

		freematrix(midT, 6 * nnod);
		freematrix(midTt, 6 * nnod);
		freematrix(midHPT, 6 * nnod);
		freematrix(midTtPtHt, 6 * nnod);

		free(lapgform);
		freematrix(lapH, 6 * nnod);
	}

	return;
}

/*
############################################################################################################################################################################################
##      ##      ####    ##          ##      ####      ######      ##                      ##          ####      ####      ##########          ####    ####          ####      ######  ######
  ######  ##  ########  ####  ####  ####  ########  ######  ######    ##  ####  ##  ####  ####  ####  ##  ######    ######  ##########  ####    ####  ######  ####    ####  ######  ########
##  ########  ######    ####  ##  ######  ########  ########  ##########  ########  ####  ####  ##  ######  ########  ################  ######  ######  ##  ######  ######  ##########  ####
##      ####            ####      ######  ########  ########      ######  ########      ######      ######      ####      ############  ######  ########  ########  ##  ##  ######      ####
########  ##  ########  ####  ##  ######  ########  ##############  ####  ########  ##  ######  ##  ############  ########  ##########  ######  ########  ########  ####    ####  ####  ####
  ######  ##  ########  ####  ##########  ####  ##  ####    ######  ####  ########  ####  ####  ########  ######    ######  ##########  ####    ########  ########  ####    ####  ######  ##
##      ##      ####    ##          ##                    ##      ####      ####      ##              ####      ####      ##########          ########      ####      ####  ##    ####
############################################################################################################################        ########################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
*/

void shellstress_DYNA(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain,
					double* iform, double* lastddisp, double* ddisp,
					double* lastudd_m, double* udd_m,
					double* lastud_m, double* ud_m,
					double* finertial, double* fdamping, double* finternal, double* fexternal,
					double alpham, double alphaf, double xi)
{
	char str[500];
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	int ngp;
	int nstress;
	int* loffset;

	double* gforminit, * eforminit;
	double* lastgform, * lasteform, * lastedisp, * lastgacc_m;
	double* gform, * eform, * edisp, * gacc_m;

	double* lasteinternal, *lasteexternal, * lastginertial_m, * lastginertial;
	double* einternal, * eexternal, * ginertial_m, * ginertial;

	double* mideinternal, * mideexternal;
	double* midginternal, * midgexternal, * midginertial;

	double** drccosinit;
	double** lastdrccos,** lastT,** lastTt,** lastHPT,** lastTtPtHt,** lastRt;
	double** drccos,** T,** Tt,** HPT,** TtPtHt,** R;
	double** midT, ** midTt, ** midHPT, ** midTtPtHt;

	double** lapH;
	double* lapgform;

	double*** C, ***B;
	double** M,** Kp;

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

		C = shellC(shell);
		B = shellB(shell);//Kp = assemshellemtx(shell);
		M = assemshellmmtx(shell);

		/*DEFORMED CONFIGURATION OF LAST LAP.*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(lastddisp, shell.node[ii]);
		}
		lastdrccos = shelldrccos(shell);
		lastgform = extractshelldisplacement(shell, lastddisp);             /*{Xg+Ug}*/
		lasteform = extractlocalcoord(lastgform,lastdrccos,nnod);           /*{Xe+Ue}*/

		lastedisp = extractdeformation(eforminit, lasteform, nnod);         /*{Ue}*/

		lastT = transmatrixIII(lastdrccos, nnod);         					/*[T]*/
		lastTt = matrixtranspose(lastT, 6 * nnod);                  		/*[Tt]*/
		lastHPT = transmatrixHPT(lasteform, lastedisp, lastT, nnod);
		lastTtPtHt = matrixtranspose(lastHPT, 6 * nnod);

		inputshell(shells, mshell, i - 1, &shell);
		lasteinternal = assemshelleinternal(&shell, B);

		lasteexternal = assemshellpvct(shell, lastdrccos);

		lastgacc_m = extractshelldisplacement(shell, lastudd_m);
		lastginertial_m = matrixvector(M, lastgacc_m, 6 * nnod);
		lastginertial = pushforward(lastgform, lastginertial_m, nnod);

		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell);
		gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
		eform = extractlocalcoord(gform,drccos,nnod); 			       	    /*{Xe+Ue}*/

		edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/

		T = transmatrixIII(drccos, nnod);
		Tt = matrixtranspose(T, 6 * nnod);
		HPT = transmatrixHPT(eform, edisp, T, nnod);
		TtPtHt = matrixtranspose(HPT, 6 * nnod);

		assemshellestrain(&shell, B, edisp);
		assemshellestress(&shell, C);
		einternal = assemshelleinternal(&shell, B);//einternal = matrixvector(Kp, edisp, 6*nnod);

		eexternal = assemshellpvct(shell, drccos);

		gacc_m = extractshelldisplacement(shell, udd_m);
		ginertial_m = matrixvector(M, gacc_m, 6 * nnod);
		ginertial = pushforward(gform, ginertial_m, nnod);

		/*(24):MID-POINT TRANSFORMATION MATRIX.*/
		midT = midpointmtx(T, lastT, alphaf, 6 * nnod);
		midTt = matrixtranspose(midT, 6 * nnod);
		midHPT = midpointmtx(HPT, lastHPT, alphaf, 6 * nnod);
		midTtPtHt = matrixtranspose(midHPT, 6 * nnod);

		/*(21)&(22)&(23):MID-POINT FORCE VECTOR.*/
		mideinternal = midpointvct(einternal, lasteinternal, alphaf-xi, 6*nnod);
		midginternal = matrixvector(midTtPtHt, mideinternal, 6 * nnod);
		midginertial = midpointvct(ginertial, lastginertial, alpham, 6*nnod);
		mideexternal = midpointvct(eexternal, lasteexternal, alphaf, 6*nnod);
		midgexternal = matrixvector(midTt, mideexternal, 6 * nnod);



		/*GLOBAL VECTOR & MATRIX ASSEMBLY.*/
		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				*(finertial + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(midginertial + 6 * ii + jj);
				*(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(midginternal + 6 * ii + jj);
				*(fexternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(midgexternal + 6 * ii + jj);
			}
		}


	//Kinetic Energy Output
	/*
	gvel_m = extractshelldisplacement(shell, ud_m);
	gmomentum_m = matrixvector(Me, gvel_m, 6 * nnod);
	for (ii = 0; ii < nnod; ii++)
	{
		for (jj = 0; jj < 3; jj++)
		{
			(mshell+i-1)->KEt += 0.5 * *(gvel_m + 6 * ii + jj) * *(gmomentum_m + 6 * ii + jj);
		}
		for (jj = 3; jj < 6; jj++)
		{
			(mshell+i-1)->KEr += 0.5 * *(gvel_m + 6 * ii + jj) * *(gmomentum_m + 6 * ii + jj);
		}
	}
	(mshell+i-1)->KE = (mshell+i-1)->KEt + (mshell+i-1)->KEr;
	free(gvel_m);
	free(gmomentum_m);
	*/



		outputshell(shells, i - 1, &shell);

		for(ii=0;ii<ngp;ii++)freematrix(*(C+ii), nstress);
		free(C);
		for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
		free(B);

		freematrix(M, 6 * nnod);

		free(loffset);

		free(eforminit);
		free(gforminit);

		free(lasteform);
		free(lastgform);
		free(lastedisp);
		free(lastgacc_m);

		free(lasteinternal);
		free(lasteexternal);
		free(lastginertial_m);
		free(lastginertial);

		free(eform);
		free(gform);
		free(edisp);
		free(gacc_m);

		free(einternal);
		free(eexternal);
		free(ginertial_m);
		free(ginertial);

		free(mideinternal);
		free(midginternal);
		free(midginertial);
		free(mideexternal);
		free(midgexternal);

		freematrix(drccosinit, 3);

		freematrix(lastdrccos, 3);
		freematrix(lastT, 6 * nnod);
		freematrix(lastTt, 6 * nnod);
		freematrix(lastHPT, 6 * nnod);
		freematrix(lastTtPtHt, 6 * nnod);

		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(TtPtHt, 6 * nnod);

		freematrix(midT, 6 * nnod);
		freematrix(midTt, 6 * nnod);
		freematrix(midHPT, 6 * nnod);
		freematrix(midTtPtHt, 6 * nnod);

	}

	return;
}









#if 0

		/*ELEMENTS ASSEMBLAGE.*/
		for (i = 1; i <= nshell; i++)
		{
			inputshell(shells, mshell, i - 1, &shell);
			nnod = shell.nnod;
			shell.sect = (shells + i - 1)->sect;
			loffset = (int*)malloc(6 * nnod * sizeof(int));
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(loffset + (6 * ii + jj)) = 6 * (shell.node[ii]->loff) + jj;
				}
			}

			/*INITIAL CONFIGURATION.*/
			for (ii = 0; ii < nnod; ii++)
			{
				inputnode(iform, shell.node[ii]);
			}
			drccosinit = shelldrccos(shell);
			gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
			eforminit = extractlocalcoord(gforminit,drccosinit,nnod);     		/*{Xe}*/

			Ke = assemshellemtx(shell);   						/*[Ke]*/
			Me = assemshellmmtx(shell);          					/*[Me]*/

			/*DEFORMED CONFIGURATION OF LAST LAP.*/
			for (ii = 0; ii < nnod; ii++)
			{
				inputnode(lastddisp, shell.node[ii]);
			}
			lastdrccos = shelldrccos(shell);
			lastgform = extractshelldisplacement(shell, lastddisp);             /*{Xg+Ug}*/
			lasteform = extractlocalcoord(lastgform,lastdrccos,nnod);           /*{Xe+Ue}*/

			lastT = transmatrixIII(lastdrccos, nnod);         					/*[T]*/
			lastTt = matrixtranspose(lastT, 6 * nnod);                  		/*[Tt]*/

			lastedisp = extractdeformation(eforminit, lasteform, nnod);         /*{Ue}*/
			lasteinternal = matrixvector(Ke, lastedisp, 6 * nnod);      		/*{Fe}=[Ke]{Ue}*/
			lastepressure = assemshellpvct(shell, lastdrccos);                	/*{Pe}*/

			lastHPT = transmatrixHPT(lasteform, lastedisp, lastT, nnod);

			lastgacc_m = extractshelldisplacement(shell, lastudd_m);
			lastginertial_m = matrixvector(Me, lastgacc_m, 6 * nnod);
			lastginertial = pushforward(lastgform, lastginertial_m, nnod);

			lastRt = pullbackmtx(lastgform, nnod);

			/*DEFORMED CONFIGURATION OF LAST ITERATION.*/
			for (ii = 0; ii < nnod; ii++)
			{
				inputnode(ddisp, shell.node[ii]);
			}
			drccos = shelldrccos(shell);
			gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
			eform = extractlocalcoord(gform,drccos,nnod);                 		/*{Xe+Ue}*/

			T = transmatrixIII(drccos, nnod);         					        /*[T]*/
			Tt = matrixtranspose(T, 6 * nnod);                  				/*[Tt]*/

			edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/
			einternal = matrixvector(Ke, edisp, 6 * nnod);      				/*{Fe}=[Ke]{Ue}*/
			epressure = assemshellpvct(shell, drccos);                		    /*{Pe}*/

			HPT = transmatrixHPT(eform, edisp, T, nnod);

			gacc_m = extractshelldisplacement(shell, udd_m);
			ginertial_m = matrixvector(Me, gacc_m, 6 * nnod);
			ginertial = pushforward(gform, ginertial_m, nnod);


			R=pushforwardmtx(gform,nnod);

			lapgform = extractshelldisplacement(shell, lapddisp);                     /*{Xg+Ug}*/
			lapH = blockjacobimtx(lapgform, NULL, NULL, nnod);

			/*(24):MID-POINT TRANSFORMATION MATRIX.*/
			midHPT = midpointmtx(HPT, lastHPT, alphaf, 6 * nnod);
			midTtPtHt = matrixtranspose(midHPT, 6 * nnod);

			midT = midpointmtx(T, lastT, alphaf, 6 * nnod);
			midTt = matrixtranspose(midT, 6 * nnod);

			/*(21)&(22)&(23):MID-POINT FORCE VECTOR.*/
			mideinternal = midpointvct(einternal, lasteinternal, alphaf-xi, 6*nnod);
			midginternal = matrixvector(midTtPtHt, mideinternal, 6 * nnod);

			midginertial = midpointvct(ginertial, lastginertial, alpham, 6*nnod);

			midepressure = midpointvct(epressure, lastepressure, alphaf, 6*nnod);
			midgpressure = matrixvector(midTt, midepressure, 6 * nnod);


			Keff=assemtmtxCR_DYNA(eform, edisp, mideinternal, T, Ke, midTtPtHt, HPT,
								  ginertial, Me, R, lastRt, lapH,
								  alphaf, alpham, xi, beta, ddt, nnod);

			symmetricmtx(Keff, 6*nnod);
			assemgstiffnessIIwithDOFelimination(gmtx, Keff, &shell, constraintmain);

			shellstress = matrixvector(DBe, edisp, 6 * nnod);
			gvel_m = extractshelldisplacement(shell, ud_m);
			gmomentum_m = matrixvector(Me, gvel_m, 6 * nnod);

			/*GLOBAL VECTOR & MATRIX ASSEMBLY.*/
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(finertial + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(midginertial + 6 * ii + jj);
					*(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(midginternal + 6 * ii + jj);
					*(fpressure + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(midgpressure + 6 * ii + jj);
					(shells+i-1)->stress[ii][jj]=*(einternal+6*ii+jj);
				}
			}




			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 3; jj++)
				{
					(mshell+i-1)->KEt += 0.5 * *(gvel_m + 6 * ii + jj) * *(gmomentum_m + 6 * ii + jj);
				}
				for (jj = 3; jj < 6; jj++)
				{
					(mshell+i-1)->KEr += 0.5 * *(gvel_m + 6 * ii + jj) * *(gmomentum_m + 6 * ii + jj);
				}
			}
			(mshell+i-1)->KE = (mshell+i-1)->KEt + (mshell+i-1)->KEr;

			free(gvel_m);
			free(gmomentum_m);
			free(shellstress);




			/*MEMORY FREE : INITIAL CONFIGURATION.*/
			free(loffset);
			freematrix(drccosinit, 3);
			free(gforminit);
			free(eforminit);
			freematrix(DBe, 6 * nnod);
			freematrix(Ke, 6 * nnod);
			freematrix(Me, 6 * nnod);

			/*MEMORY FREE : DEFORMED CONFIGURATION OF LAST LAP.*/
			freematrix(lastdrccos, 3);
			freematrix(lastT, 6 * nnod);
			freematrix(lastTt, 6 * nnod);
			free(lastgform);
			free(lasteform);
			free(lastedisp);
			free(lasteinternal);
			freematrix(lastHPT, 6 * nnod);
			free(lastgacc_m);
			freematrix(lastRt, 6 * nnod);
			free(lastginertial_m);
			free(lastginertial);

			free(lastepressure);

			/*MEMORY FREE : DEFORMED CONFIGURATION OF LAST ITERATION.*/
			freematrix(drccos, 3);
			freematrix(T, 6 * nnod);
			freematrix(Tt, 6 * nnod);
			free(gform);
			free(eform);
			free(edisp);
			free(einternal);
			freematrix(HPT, 6 * nnod);
			free(gacc_m);
			freematrix(R, 6 * nnod);
			free(ginertial_m);
			free(ginertial);

			free(lapgform);
			freematrix(lapH, 6 * nnod);

			free(epressure);

			/*MEMORY FREE : MID-POINT.*/

			freematrix(midHPT, 6 * nnod);
			freematrix(midTtPtHt, 6 * nnod);
			freematrix(midT, 6 * nnod);
			freematrix(midTt, 6 * nnod);
			free(mideinternal);
			free(midginternal);
			free(midginertial);
			free(midepressure);
			free(midgpressure);
			freematrix(Keff, 6 * nnod);

		}



#endif



