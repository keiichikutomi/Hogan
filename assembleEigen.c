
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
void assemelemEigen(struct owire* elems, struct memoryelem* melem, int nelem, long int* constraintmain,struct oconf* confs,
					double* dirichletdisp, double* dirichletvct,
					std::vector<Triplet>& Mtriplet,std::vector<Triplet>& Ktriplet,
					double* iform, double* ddisp,
					int SYMFLAG)
{
	struct owire elem;
	int i,j,ii,jj;
	int nnod;
	int* loffset;

	double* gforminit, * eforminit;
	double* gform, * eform, * edisp;

	double* ginternal, * einternal;

	double** drccosinit;
	double** drccos, ** T, **HPT;

	double** Kp,** Kt;



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

		Kp = assememtx(elem);
		//Kp = modifyhinge(elem,Kp);             /*MODIFY MATRIX.*/

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
		HPT = transmatrixHPT(eform, edisp, T, nnod);

		einternal = matrixvector(Kp, edisp, 6 * nnod);          			/*{Fe}=[Ke]{Ue}.*/

		//Kt = assemtmtxCR(Ke, eform, edisp, einternal, T, HPT, nnod);	/*TANGENTIAL MATRIX[Kt].*/
		Kp = transformationIII(Kp, HPT, 6*nnod);/*[Ke]=[Tt][Pt][Ht][K][H][P][T]*/
		Kt = assemgmtxCR(eform, edisp, einternal, T, nnod);/*[Kg]=[Kgr]+[Kgp]+[Kgm]*/
		for (ii = 0; ii < 6*nnod; ii++)
		{
			for (jj = 0; jj < 6*nnod; jj++)
			{
				*(*(Kt + ii) + jj) += *(*(Kp + ii) + jj);/*[Kt]=[Ke]+[Kg]*/
			}
		}


		if(SYMFLAG==1)
		{
		  symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
		}

		assemtripletelem(Ktriplet, Kt, &elem, constraintmain, confs);


		/*DIRICHLET*/
		double* gdirichlet,* ginternal;
		if(dirichletdisp!=NULL && dirichletvct!=NULL)
		{
			gdirichlet = extractdisplacement(elem, dirichletdisp);
			ginternal = matrixvector(Kt,gdirichlet,6*nnod);
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(dirichletvct + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				}
			}
			free(gdirichlet);
			free(ginternal);
		}



		free(loffset);
		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);

		free(einternal);

		freematrix(drccosinit, 3);
		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(Kp, 6 * nnod);
		freematrix(Kt, 6 * nnod);


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

void assemshellEigen(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain,struct oconf* confs,
					 double* dirichletdisp, double* dirichletvct,
					 std::vector<Triplet>& Mtriplet,std::vector<Triplet>& Ktriplet,
					 double* iform, double* ddisp,
					 int SYMFLAG)
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
		C = shellCconsistent(shell);
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
		if(SYMFLAG==1)
		{
		  symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
		}

		assemtripletshell(Ktriplet, Kt, &shell, constraintmain, confs);

		M = transformationIII(M, T, 6*nnod);
		assemtripletshell(Mtriplet, M, &shell, constraintmain, confs);

		/*DIRICHLET*/
		double* gdirichlet,* ginternal;
		if(dirichletdisp!=NULL && dirichletvct!=NULL)
		{
			gdirichlet = extractshelldisplacement(shell, dirichletdisp);
			ginternal = matrixvector(Kt,gdirichlet,6*nnod);
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(dirichletvct + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				}
			}
			free(gdirichlet);
			free(ginternal);
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





void assemelemEigen_DYNA(struct owire* elems, struct memoryelem* melem, int nelem, long int* constraintmain,struct oconf* confs,
						 double* dirichletdisp, double* dirichletvct,
						 std::vector<Triplet>& Ktriplet,std::vector<Triplet>& Ktriplet2,
						 double* iform, double* lastddisp, double* ddisp,
						 double* ud_m, double* udd_m,
						 double alpham, double alphaf, double xi, double beta, double ddt, int SYMFLAG)
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

		assemtripletelem(Ktriplet2, Kint, &elem, constraintmain, confs);

		Keff = assemtmtxCR_DYNA(Kint,
								ginertial, M, R, lastRt, lapH,
								alpham, beta, ddt, nnod);
		if(SYMFLAG==1)
		{
		  symmetricmtx(Keff, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
		}
		assemtripletelem(Ktriplet, Keff, &elem, constraintmain, confs);

		/*DIRICHLET*/
		double* gdirichlet,* ginternal;
		if(dirichletdisp!=NULL && dirichletvct!=NULL)
		{
			gdirichlet = extractdisplacement(elem, dirichletdisp);
			ginternal = matrixvector(Keff,gdirichlet,6*nnod);
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(dirichletvct + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				}
			}
			free(gdirichlet);
			free(ginternal);
		}


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
void assemshellEigen_DYNA(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain,struct oconf* confs,
						  double* dirichletdisp, double* dirichletvct,
						  std::vector<Triplet>& Ktriplet,std::vector<Triplet>& Ktriplet2,
						  double* iform, double* lastddisp, double* ddisp,
						  double* ud_m, double* udd_m,
						  double alpham, double alphaf, double xi, double beta, double ddt, int SYMFLAG)
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
		C = shellCconsistent(shell);
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

		assemtripletshell(Ktriplet2, Kint, &shell, constraintmain, confs);


		Keff = assemtmtxCR_DYNA(Kint,
								ginertial, M, R, lastRt, lapH,
								alpham, beta, ddt, nnod);
		if(SYMFLAG==1)
		{
		  symmetricmtx(Keff, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
		}
		assemtripletshell(Ktriplet, Keff, &shell, constraintmain, confs);

		/*DIRICHLET*/
		double* gdirichlet,* ginternal;
		if(dirichletdisp!=NULL && dirichletvct!=NULL)
		{
			gdirichlet = extractshelldisplacement(shell, dirichletdisp);
			ginternal = matrixvector(Keff,gdirichlet,6*nnod);
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(dirichletvct + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				}
			}
			free(gdirichlet);
			free(ginternal);
		}

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



