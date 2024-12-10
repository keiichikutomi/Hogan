
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

		Kt = assemtmtxCR(Ke, eform, edisp, einternal, T, HPT, nnod);	/*TANGENTIAL MATRIX[Kt].*/
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
	double** drccos,** T,** HPT;

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

		B = shellB(shell);
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
		HPT = transmatrixHPT(eform, edisp, T, nnod);

		inputshell(shells, NULL, i - 1, &shell);
		C = shellCconsistentilyushin(shell);
		einternal = assemshelleinternal(&shell, B);

		Kp = assemshellpmtx(shell,C,B);
		Kp = transformationIII(Kp, HPT, 6*nnod);/*[Ke]=[Tt][Pt][Ht][K][H][P][T]*/
		Kt = assemgmtxCR(eform, edisp, einternal, T, nnod);/*[Kg]=[Kgr]+[Kgp]+[Kgm]*/
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
		freematrix(HPT, 6 * nnod);
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
		ginternal = matrixvector(TtPtHt, einternal, 6 * nnod);

		eexternal = assemshellpvct(shell, drccos);
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




void assemelem_DYNA(struct owire* elems, struct memoryelem* melem, int nelem, long int* constraintmain,
					 struct gcomponent* gmtx,struct gcomponent* gmtx2,
					 double* iform, double* lastddisp, double* ddisp,
					 double* ud_m, double* udd_m,
					 double alpham, double alphaf, double xi, double beta, double ddt)
{
	char str[500];
	struct owire elem;
	int i,j,ii,jj;
	int nnod;
	int* loffset;

	double* gforminit, * eforminit;
	double* lastgform, * lasteform, * lastedisp, * lastgacc_m;
	double* gform, * eform, * edisp, * gacc_m;

	double* lasteinternal;
	double* einternal, * ginertial_m, * ginertial;
	double* mideinternal;

	double** drccosinit;
	double** lastdrccos,** lastT,** lastHPT,** lastRt;
	double** drccos,** T,** HPT,** R;
	double** midHPT, ** midTtPtHt;

	double** lapH;
	double* lapgform;
	double** M, ** Kp, ** Kint, ** Keff;

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


		//K = assememtx(elem);
		//M = assemmmtx(elem, drccosinit);          					/*[Me]*/

		/*DEFORMED CONFIGURATION OF LAST LAP.*/
		/*DEFORMED CONFIGURATION OF LAST LAP.*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(lastddisp, elem.node[ii]);
		}
		lastgform = extractdisplacement(elem, lastddisp);
		lastdrccos = updatedrccos(drccosinit, gforminit, lastgform);
		lasteform = extractlocalcoord(lastgform,lastdrccos,nnod);

		lastedisp = extractdeformation(eforminit, lasteform, nnod);         /*{Ue}*/

		lastT = transmatrixIII(lastdrccos, nnod);         					/*[T]*/
		lastHPT = transmatrixHPT(lasteform, lastedisp, lastT, nnod);//
		lastRt = pullbackmtx(lastgform, nnod);

		inputelem(elems,melem,i-1,&elem);
		//lasteinternal = assemeinternal(&elem);

		/*DEFORMED CONFIDURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, elem.node[ii]);
		}
		gform = extractdisplacement(elem, ddisp);
		drccos = updatedrccos(drccosinit, gforminit, gform);
		eform = extractlocalcoord(gform,drccos,nnod);

		edisp = extractdeformation(eforminit, eform, nnod);//

		T = transmatrixIII(drccos, nnod);//         							/*[T].*/
		HPT = transmatrixHPT(eform, edisp, T, nnod);//
		R = pushforwardmtx(gform,nnod);

		inputelem(elems,NULL,i-1,&elem); /*CHECK FUNCTION*/
		//C = shellCconsistentilyushin(shell);
		//einternal = assemeinternal(&elem);//

		gacc_m = extractdisplacement(elem, udd_m);
		ginertial_m = matrixvector(M, gacc_m, 6 * nnod);
		ginertial = pushforward(gform, ginertial_m, nnod);//

		//Kp = assempmtx(shell);

		lapgform = extractdeformation(lastgform, gform, nnod);
		lapH = blockjacobimtx(lapgform, NULL, NULL, nnod);



		/*(24):MID-POINT TRANSFORMATION MATRIX.*/
		midHPT = midpointmtx(HPT, lastHPT, alphaf, 6 * nnod);
		midTtPtHt = matrixtranspose(midHPT, 6 * nnod);

		/*(21)&(22)&(23):MID-POINT FORCE VECTOR.*/
		mideinternal = midpointvct(einternal, lasteinternal, alphaf-xi, 6*nnod);

		Kint = assemtmtxCR_MID(Kp, eform, edisp, mideinternal, T, midTtPtHt, HPT, alphaf, xi, nnod);
		symmetricmtx(Kint, 6*nnod);
		if(gmtx2!=NULL)assemgstiffnesswithDOFelimination(gmtx2, Kint, &elem, constraintmain);

		Keff = assemtmtxCR_DYNA(Kint,
								ginertial, M, R, lastRt, lapH,
								alpham, beta, ddt, nnod);
		symmetricmtx(Keff, 6*nnod);
		assemgstiffnesswithDOFelimination(gmtx, Keff, &elem, constraintmain);

		//freematrix(M, 6 * nnod);
		//freematrix(Kp, 6 * nnod);
		freematrix(Kint, 6 * nnod);
		freematrix(Keff, 6 * nnod);

		free(loffset);

		free(eforminit);
		free(gforminit);

		free(lasteform);
		free(lastgform);
		free(lastedisp);
		free(lasteinternal);

		free(eform);
		free(gform);
		free(edisp);
		free(einternal);
		free(gacc_m);
		free(ginertial_m);
		free(ginertial);

		free(mideinternal);



		freematrix(drccosinit, 3);

		freematrix(lastdrccos, 3);
		freematrix(lastT, 6 * nnod);
		freematrix(lastHPT, 6 * nnod);
		freematrix(lastRt, 6 * nnod);

		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(R, 6 * nnod);

		freematrix(midHPT, 6 * nnod);
		freematrix(midTtPtHt, 6 * nnod);

		free(lapgform);
		freematrix(lapH, 6 * nnod);
	}

	return;
}




void elemstress_DYNA(struct owire* elems, struct memoryelem* melem, int nelem, long int* constraintmain,
					 double* iform, double* lastddisp, double* ddisp,
					 double* lastudd_m, double* udd_m,
					 double* lastud_m, double* ud_m,
					 double* finertial, double* fdamping, double* finternal, double* fexternal,
					 double alpham, double alphaf, double xi)
{
	char str[500];
	struct owire elem;
	int i,j,ii,jj;
	int nnod;
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

	double** M;


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

		//K = assememtx(elem);
		//M = assemmmtx(elem, drccosinit);          					/*[Me]*/

		/*DEFORMED CONFIGURATION OF LAST LAP.*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(lastddisp, elem.node[ii]);
		}
		lastgform = extractdisplacement(elem, lastddisp);
		lastdrccos = updatedrccos(drccosinit, gforminit, lastgform);
		lasteform = extractlocalcoord(lastgform,lastdrccos,nnod);

		lastedisp = extractdeformation(eforminit, lasteform, nnod);

		lastT = transmatrixIII(lastdrccos, nnod);
		lastTt = matrixtranspose(lastT, 6 * nnod);
		lastHPT = transmatrixHPT(lasteform, lastedisp, lastT, nnod);
		lastTtPtHt = matrixtranspose(lastHPT, 6 * nnod);

		//lasteinternal = assemeinternal(&elem);
		//lasteexternal = assempvct(elem, lastdrccos);

		lastgacc_m = extractdisplacement(elem, lastudd_m);
		lastginertial_m = matrixvector(M, lastgacc_m, 6 * nnod);
		lastginertial = pushforward(lastgform, lastginertial_m, nnod);

		/*DEFORMED CONFIDURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, elem.node[ii]);
		}
		gform = extractdisplacement(elem, ddisp);
		drccos = updatedrccos(drccosinit, gforminit, gform);
		eform = extractlocalcoord(gform,drccos,nnod);

		edisp = extractdeformation(eforminit, eform, nnod);

		T = transmatrixIII(drccos, nnod);
		Tt = matrixtranspose(T, 6 * nnod);
		HPT = transmatrixHPT(eform, edisp, T, nnod);
		TtPtHt = matrixtranspose(HPT, 6 * nnod);

		//assemshellestrain(&shell, B, edisp);
		//assemshellestress(&shell, C);
		//einternal = assemeinternal(&elem);

		//eexternal = assempvct(elem, drccos);

		gacc_m = extractdisplacement(elem, udd_m);
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





		//outputelem(elems, i - 1, &elem);


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
					 struct gcomponent* gmtx,struct gcomponent* gmtx2,
					 double* iform, double* lastddisp, double* ddisp,
					 double* ud_m, double* udd_m,
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

	double* lasteinternal;
	double* einternal, * ginertial_m, * ginertial;
	double* mideinternal;

	double** drccosinit;
	double** lastdrccos,** lastT,** lastHPT,** lastRt;
	double** drccos,** T,** HPT,** R;
	double** midHPT, ** midTtPtHt;

	double** lapH;
	double* lapgform;
	double*** C, ***B;
	double** M, ** Kp, ** Kint, ** Keff;

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
		drccosinit = shelldrccos(shell);//
		gforminit = extractshelldisplacement(shell, iform);//
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);//

		B = shellB(shell);//
		M = assemshellmmtx(shell);         				                   /*[Me]*/

		/*DEFORMED CONFIGURATION OF LAST LAP.*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(lastddisp, shell.node[ii]);
		}
		lastdrccos = shelldrccos(shell); //
		lastgform = extractshelldisplacement(shell, lastddisp);//             /*{Xg+Ug}*/
		lasteform = extractlocalcoord(lastgform,lastdrccos,nnod); //          /*{Xe+Ue}*/

		lastedisp = extractdeformation(eforminit, lasteform, nnod);           		/*{Ue}*/

		lastT = transmatrixIII(lastdrccos, nnod);//         					/*[T]*/
		lastHPT = transmatrixHPT(lasteform, lastedisp, lastT, nnod);//
		lastRt = pullbackmtx(lastgform, nnod);

		inputshell(shells, mshell, i - 1, &shell);
		lasteinternal = assemshelleinternal(&shell, B);//


		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell);//
		gform = extractshelldisplacement(shell, ddisp);//
		eform = extractlocalcoord(gform,drccos,nnod);// 			       	    /*{Xe+Ue}*/

		edisp = extractdeformation(eforminit, eform, nnod);//

		T = transmatrixIII(drccos, nnod);//         							/*[T].*/
		HPT = transmatrixHPT(eform, edisp, T, nnod);//
		R = pushforwardmtx(gform,nnod);

		inputshell(shells, NULL, i - 1, &shell);
		C = shellCconsistentilyushin(shell);
		einternal = assemshelleinternal(&shell, B);//

		gacc_m = extractshelldisplacement(shell, udd_m);
		ginertial_m = matrixvector(M, gacc_m, 6 * nnod);
		ginertial = pushforward(gform, ginertial_m, nnod);//

		Kp = assemshellpmtx(shell,C,B);

		lapgform = extractdeformation(lastgform, gform, nnod);
		lapH = blockjacobimtx(lapgform, NULL, NULL, nnod);

		/*(24):MID-POINT TRANSFORMATION MATRIX.*/
		midHPT = midpointmtx(HPT, lastHPT, alphaf, 6 * nnod);
		midTtPtHt = matrixtranspose(midHPT, 6 * nnod);

		/*(21)&(22)&(23):MID-POINT FORCE VECTOR.*/
		mideinternal = midpointvct(einternal, lasteinternal, alphaf-xi, 6*nnod);//

		Kint = assemtmtxCR_MID(Kp, eform, edisp, mideinternal, T, midTtPtHt, HPT, alphaf, xi, nnod);
		symmetricmtx(Kint, 6*nnod);
		if(gmtx2!=NULL)assemgstiffnessIIwithDOFelimination(gmtx2, Kint, &shell, constraintmain);

		Keff = assemtmtxCR_DYNA(Kint,
								ginertial, M, R, lastRt, lapH,
								alpham, beta, ddt, nnod);
		symmetricmtx(Keff, 6*nnod);
		assemgstiffnessIIwithDOFelimination(gmtx, Keff, &shell, constraintmain);

		for(ii=0;ii<ngp;ii++)freematrix(*(C+ii), nstress);
		free(C);
		for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
		free(B);

		freematrix(M, 6 * nnod);
		freematrix(Kp, 6 * nnod);
		freematrix(Kint, 6 * nnod);
		freematrix(Keff, 6 * nnod);

		free(loffset);

		free(eforminit);
		free(gforminit);

		free(lasteform);
		free(lastgform);
		free(lastedisp);
		free(lasteinternal);

		free(eform);
		free(gform);
		free(edisp);
		free(einternal);
		free(gacc_m);
		free(ginertial_m);
		free(ginertial);

		free(mideinternal);



		freematrix(drccosinit, 3);

		freematrix(lastdrccos, 3);
		freematrix(lastT, 6 * nnod);
		freematrix(lastHPT, 6 * nnod);
		freematrix(lastRt, 6 * nnod);

		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(R, 6 * nnod);

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
	double** M;

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
		einternal = assemshelleinternal(&shell, B);

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


void strainenergy(struct arclmframe* af, double* Wet, double* Wpt)
{
	int i,j;

	*Wet = 0.0;
	*Wpt = 0.0;

	for (i = 0; i < af->nelem; i++)
	{
		for (j = 0; j < 2; j++)
		{
			*Wet += (af->elems+i)->Ee[j];
			*Wpt += (af->elems+i)->Ep[j];
		}
	}
	for (i = 0; i < af->nshell; i++)
	{
		for (j = 0; j < (af->shells+i)->ngp; j++)
		{
			*Wet += ((af->shells+i)->gp[j]).Ee;
			*Wpt += ((af->shells+i)->gp[j]).Ep;
		}
	}

	return;
}



void kineticenergy(struct arclmframe* af, double* ud_m, double* Wkt)
{
	int i,j,k;
	int nnod,ngp,nstress;
	double** M;
	double* gvel_m, * gmomentum_m;

	*Wkt = 0.0;

	for(i = 0; i < af->nshell; i++)
	{
		nnod = (af->shells+i)->nnod;
		M = assemshellmmtx(*(af->shells+i));
		gvel_m = extractshelldisplacement(*(af->shells+i), ud_m);
		gmomentum_m = matrixvector(M, gvel_m, 6 * nnod);
		for (j = 0; j < nnod; j++)
		{
			for (k= 0; k < 6; k++)
			{
				*Wkt += 0.5 * *(gvel_m + 6 * j + k) * *(gmomentum_m + 6 * j + k);
			}
		}
		free(gvel_m);
		free(gmomentum_m);
		freematrix(M,6*nnod);
	}

	for (i = 0; i < af->nnode; i++)
	{
		for (j = 0; j < 3; j++)
		{
			*Wkt += 0.5 * (*(af->nmass + i)) * (*(ud_m + 6 * i + j)) *(*(ud_m + 6 * i + j));
			//Wot -= 0.5 * (*(nmass + ii)) * (2 * tacc[jj] - dacc[jj]) *(*(du + 6 * ii + jj));
		}
	}

	return;
}















