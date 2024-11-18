/*FOR ROTATION CALCURATION*/
double* rotationvct(double** rmtx);
double** rotationmtx(double* rvct);/*[R]*/
double** spinmtx(double* rvct);
double** jacobimtx(double* rvct);/*[H]*/

double** spinfittermtx(double* eform, int nnod);/*[G]*/
double** projectionmtx(double* eform, double** G, int nnod);/*[P]*/
double** blockjacobimtx(double* edisp, double* estress, double** M, int nnod);/*[H]&[M]*/

double** transmatrixHPT(double* eform, double* edisp, double** T, int nnod);
double** assemtmtxCR(double** Ke, double* eform, double* edisp, double* estress, double* gstress, double** T, double** HPT, int nnod);
double** assemgmtxCR(double* eform, double* edisp, double* estress, double* gstress, double** T, double** HPT, int nnod);
void symmetricmtx(double** estiff, int msize);

void updaterotation(double* ddisp, double* gvct, int nnode);
double* quaternionvct(double* rvct);
double** updatedrccos(double** drccosinit, double* gforminit, double* gform);
double** interpolatermtx(double* rvct1, double* rvct2, double alpha);

double* extractlocalcoord(double* gform, double** drccos, double nnod);
double* extractdeformation(double* eforminit, double* eform, int nnod);

double* pullback(double* ddisp, double* gvct_s, int nnode);
double** pullbackmtx(double* gform, int nnod);
double* pushforward(double* ddisp, double* gvct_m, int nnode);
double** pushforwardmtx(double* gform, int nnod);
double* midpointvct(double* vct,double* lastvct,double alpha,int size);
double** midpointmtx(double** mtx,double** lastmtx,double alpha,int size);






double* rotationvct(double** rmtx)
{
	char str[500];
	double c, s;
	double* rvct;
	double theta;
	rvct = (double*)malloc(3 * sizeof(double));
	c = 0.5*(  *(*(rmtx + 0) + 0) + *(*(rmtx + 1) + 1) + *(*(rmtx + 2) + 2) - 1  );/*cos(theta)*/
	s = 0.5*(sqrt(pow( *(*(rmtx + 2) + 1) - *(*(rmtx + 1) + 2) ,2)
				+ pow( *(*(rmtx + 0) + 2) - *(*(rmtx + 2) + 0) ,2)
				+ pow( *(*(rmtx + 1) + 0) - *(*(rmtx + 0) + 1) ,2))); /*|sin(theta)|>=0*/

	theta = atan2(s , c);/*theta>=0*/

	if(s > 1e-8)
	{
		*(rvct + 0) = 0.5 * ((*(*(rmtx + 2) + 1) - *(*(rmtx + 1) + 2)))* theta / s;
		*(rvct + 1) = 0.5 * ((*(*(rmtx + 0) + 2) - *(*(rmtx + 2) + 0)))* theta / s;
		*(rvct + 2) = 0.5 * ((*(*(rmtx + 1) + 0) - *(*(rmtx + 0) + 1)))* theta / s;
	}
	else if(c < 0.0)/*theta=PI*/
	{
		if (*(*(rmtx + 0) + 0) >= *(*(rmtx + 1) + 1) && *(*(rmtx + 0) + 0) >= *(*(rmtx + 2) + 2))
		{
			*(rvct + 0) = sqrt( ( *(*(rmtx + 0) + 0) + 1.0 ) / 2.0 );
			*(rvct + 1) = 0.5 * *(*(rmtx + 0) + 1) / (*(rvct + 0));
			*(rvct + 2) = 0.5 * *(*(rmtx + 0) + 2) / (*(rvct + 0));
		}
		if (*(*(rmtx + 1) + 1) >= *(*(rmtx + 0) + 0) && *(*(rmtx + 1) + 1) >= *(*(rmtx + 2) + 2))
		{
			*(rvct + 1) = sqrt( ( *(*(rmtx + 1) + 1) + 1.0 ) / 2.0 );
			*(rvct + 0) = 0.5 * *(*(rmtx + 1) + 0) / (*(rvct + 1));
			*(rvct + 2) = 0.5 * *(*(rmtx + 1) + 2) / (*(rvct + 1));
		}
		if (*(*(rmtx + 2) + 2) >= *(*(rmtx + 0) + 0) && *(*(rmtx + 2) + 2) >= *(*(rmtx + 1) + 1))
		{
			*(rvct + 2) = sqrt( ( *(*(rmtx + 2) + 2) + 1.0 ) / 2.0 );
			*(rvct + 0) = 0.5 * *(*(rmtx + 2) + 0) / (*(rvct + 2));
			*(rvct + 1) = 0.5 * *(*(rmtx + 2) + 1) / (*(rvct + 2));
		}

		*(rvct + 0) *= PI;
		*(rvct + 1) *= PI;
		*(rvct + 2) *= PI;

		//sprintf(str,"%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",s,c,theta,*(rvct + 0),*(rvct + 1),*(rvct + 2), *(*(rmtx + 0) + 0), *(*(rmtx + 1) + 1), *(*(rmtx + 2) + 2));
		//errormessage(str);

	}
	else/*theta=0*/
	{
		*(rvct + 0) = 0.0;
		*(rvct + 1) = 0.0;
		*(rvct + 2) = 0.0;
	}

	return rvct;
}


double** rotationmtx(double* rvct)
{
	int i;
	double** rmtx;
	double* n;
	double theta;

	n = (double*)malloc(3 * sizeof(double));
	rmtx = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(rmtx + i) = (double*)malloc(3 * sizeof(double));
	}
	theta = sqrt(*(rvct + 0) * *(rvct + 0) + *(rvct + 1) * *(rvct + 1) + *(rvct + 2) * *(rvct + 2));
	if (theta > 1e-15)
	{
		for (i = 0; i < 3; i++)
		{
			*(n + i) = *(rvct + i) / theta;
		}

		*(*(rmtx + 0) + 0) =  cos(theta)            + (1 - cos(theta)) * *(n + 0) * *(n + 0);
		*(*(rmtx + 0) + 1) = -sin(theta) * *(n + 2) + (1 - cos(theta)) * *(n + 0) * *(n + 1);
		*(*(rmtx + 0) + 2) =  sin(theta) * *(n + 1) + (1 - cos(theta)) * *(n + 0) * *(n + 2);
		*(*(rmtx + 1) + 0) =  sin(theta) * *(n + 2) + (1 - cos(theta)) * *(n + 1) * *(n + 0);
		*(*(rmtx + 1) + 1) =  cos(theta)            + (1 - cos(theta)) * *(n + 1) * *(n + 1);
		*(*(rmtx + 1) + 2) = -sin(theta) * *(n + 0) + (1 - cos(theta)) * *(n + 1) * *(n + 2);
		*(*(rmtx + 2) + 0) = -sin(theta) * *(n + 1) + (1 - cos(theta)) * *(n + 2) * *(n + 0);
		*(*(rmtx + 2) + 1) =  sin(theta) * *(n + 0) + (1 - cos(theta)) * *(n + 2) * *(n + 1);
		*(*(rmtx + 2) + 2) =  cos(theta)            + (1 - cos(theta)) * *(n + 2) * *(n + 2);
	}
	else
	{
		*(*(rmtx + 0) + 0) = 1.0;
		*(*(rmtx + 0) + 1) = 0.0;
		*(*(rmtx + 0) + 2) = 0.0;
		*(*(rmtx + 1) + 0) = 0.0;
		*(*(rmtx + 1) + 1) = 1.0;
		*(*(rmtx + 1) + 2) = 0.0;
		*(*(rmtx + 2) + 0) = 0.0;
		*(*(rmtx + 2) + 1) = 0.0;
		*(*(rmtx + 2) + 2) = 1.0;
	}
	free(n);
	return rmtx;
}

double** spinmtx(double* rvct)
{
	int i;
	double** smtx;

	smtx = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(smtx + i) = (double*)malloc(3 * sizeof(double));
	}

	*(*(smtx + 0) + 0) = 0.0;
	*(*(smtx + 0) + 1) = -*(rvct + 2);
	*(*(smtx + 0) + 2) = *(rvct + 1);
	*(*(smtx + 1) + 0) = *(rvct + 2);
	*(*(smtx + 1) + 1) = 0.0;
	*(*(smtx + 1) + 2) = -*(rvct + 0);
	*(*(smtx + 2) + 0) = -*(rvct + 1);
	*(*(smtx + 2) + 1) = *(rvct + 0);
	*(*(smtx + 2) + 2) = 0.0;

	return smtx;
}

double** spinfittermtx(double* eform, int nnod)
/*spin-fitter matrix for beam & shell*/
/*G     :3*6nnod matrix(rotational variation by translational motion)*/
/*input :6nnod vector(variation of nodal displacement in local which inclouds noneffective rotational DOF)*/
/*output:3 vector(variation of local coord psuedo rotation)*/
{
	int i, j;
	double A, len;
	double x1, x2, y1, y2, z1;
	double** G;


	G = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(G + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(G + i) + j) = 0.0;
		}
	}


	if(nnod==2)/*[G1,0,G2,0]*/
	{

		x1 = *(eform + 6) - *(eform + 0);
		y1 = *(eform + 7) - *(eform + 1);
		z1 = *(eform + 8) - *(eform + 2);
		len = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

		/*G1*/
		*(*(G + 1) + 2) =  1 / len;
		*(*(G + 2) + 1) = -1 / len;
		/*G2*/
		*(*(G + 1) + 8) = -1 / len;
		*(*(G + 2) + 7) =  1 / len;
	}
	if(nnod==3)/*[G1,0,G2,0,G3,0]*/
	{

		x1 = *(eform + 6) - *(eform + 0);
		y1 = *(eform + 7) - *(eform + 1);
		x2 = *(eform + 12) - *(eform + 0);
		y2 = *(eform + 13) - *(eform + 1);
		A = 0.5 * (x1 * y2 - x2 * y1);
		len = sqrt(x1 * x1 + y1 * y1);

		/*G1*/
		*(*(G + 0) + 2) = 0.5 * (*(eform + 12) - *(eform + 6)) / A;
		*(*(G + 1) + 2) = 0.5 * (*(eform + 13) - *(eform + 7)) / A;
		*(*(G + 2) + 1) = -1 / len;
		/*G2*/
		*(*(G + 0) + 8) = 0.5 * (*(eform + 0) - *(eform + 12)) / A;
		*(*(G + 1) + 8) = 0.5 * (*(eform + 1) - *(eform + 13)) / A;
		*(*(G + 2) + 7) = 1 / len;
		/*G3*/
		*(*(G + 0) + 14) = 0.5 * (*(eform + 6) - *(eform + 0)) / A;
		*(*(G + 1) + 14) = 0.5 * (*(eform + 7) - *(eform + 1)) / A;
	}
	return G;
}

double** projectionmtx(double* eform, double** G,int nnod)
/*element projection matrix for beam & shell*/
/*G     :6nnod*6nnod matrix(variation correction)*/
/*input :6nnod vector(initial variation of displacement in local coord)*/
/*output:6nnod vector(corrected variation of displacement in local coord)*/
{
	int i, j, a, b;
	double* node;
	double** P, ** Gu, ** S, ** SGu;

	P = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(P + i) = (double*)malloc(6*nnod * sizeof(double));
	}

	/*eform : latest coordination of nodes in local*/
	node = (double*)malloc(3 * sizeof(double));

	/*for extracting Gu_b(b=1,2,3) from G*/
	Gu = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(Gu + i) = (double*)malloc(3 * sizeof(double));
	}
	//Pab = [ delta*I-1/3*I+S_a*Gu_b    0       ]
	//      [         -Gu_b             delta*I ]

	for (a = 0; a < nnod; a++)/*row 6*a+0~6*a+5 of P*/
	{
		for (i = 0; i < 3; i++)
		{
			*(node + i) = *(eform + 6 * a + i);
		}
		S = spinmtx(node);
		for (b = 0; b < nnod; b++)/*column 6*b+0~6*b+5 of P*/
		{
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(Gu + i) + j) = *(*(G + i) + 6 * b + j);
				}
			}
			SGu = matrixmatrix(S, Gu, 3);
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(P + 6 * a + 3 + i) + 6 * b + j) = -*(*(Gu + i) + j);/*-Gu_b*/
					*(*(P + 6 * a + i) + 6 * b + 3 + j) = 0.0;				/*0*/
					if (a == b && i == j)
					{
						*(*(P + 6 * a + i) + 6 * b + j) = 1.0 - 1.0/(double)nnod;	/*delta*I-1/3*I*/
						*(*(P + 6 * a + 3 + i) + 6 * b + 3 + j) = 1.0;		/*delta*I*/
					}
					else if (i == j)
					{
						*(*(P + 6 * a + i) + 6 * b + j) = - 1.0/(double)nnod;		/*delta*I-1/3*I*/
						*(*(P + 6 * a + 3 + i) + 6 * b + 3 + j) = 0.0;		/*delta*I*/
					}
					else
					{
						*(*(P + 6 * a + i) + 6 * b + j) = 0.0;				/*delta*I-1/3*I*/
						*(*(P + 6 * a + 3 + i) + 6 * b + 3 + j) = 0.0;		/*delta*I*/
					}
					*(*(P + 6 * a + i) + 6 * b + j) += *(*(SGu + i) + j);	/*S_a*Gu_b*/
				}
			}
			freematrix(SGu, 3);
		}
		freematrix(S, 3);
	}
	free(node);
	freematrix(Gu, 3);
	return P;
}

double** jacobimtx(double* rvct)
/*transformation from additional infinitesimal incremental rotation to increment of total rotational pseudo-vector*/
{
	int n, i, j;
	double theta, eta;
	double** thetaspin, ** thetaspin2;
	double** H;


	H = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(H + i) = (double*)malloc(3 * sizeof(double));
		for (j = 0; j < 3; j++)
		{
			*(*(H + i) + j) = 0.0;
		}
	}

	theta = vectorlength(rvct,3);

	if (theta < PI / 30.0)
	{
		eta = 1.0 / 12.0 + pow(theta, 2) / 720.0 + pow(theta, 4) / 30240.0 + pow(theta, 6) / 1209600.0;
	}
	else
	{
		eta = (1.0 - 0.5 * theta * (1.0 / tan(0.5 * theta))) / pow(theta, 2);
	}
	thetaspin = spinmtx(rvct);
	thetaspin2 = matrixmatrix(thetaspin, thetaspin, 3);


	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			*(*(H + i) + j) = -0.5 * *(*(thetaspin + i) + j) + eta * *(*(thetaspin2 + i) + j);
			if (i == j)*(*(H + i) + j) += 1.0;
		}
	}
	freematrix(thetaspin, 3);
	freematrix(thetaspin2, 3);
	return H;
}

double** blockjacobimtx(double* edisp, double* estress, double** M, int nnod)
{
	int n, i, j;
	double theta, dot, eta, mu;
	double* rvct, * mvct;
	double** thetaspin, ** thetaspin2, ** mspin, ** mtheta, ** mtheta2;
	double** H, ** Ha, ** Ma;

	/*H=diag[I, H_1, I, H_2(, I, H_3)]*/
	H = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(H + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(H + i) + j) = 0.0;
		}
	}
	rvct = (double*)malloc(3 * sizeof(double));

	if(estress!=NULL && M!=NULL)
	{
		mtheta = (double**)malloc(3 * sizeof(double*));                         /*{ma}{θat}*/
		for (i = 0; i < 3; i++)
		{
			*(mtheta + i) = (double*)malloc(3 * sizeof(double));
		}
		mvct = (double*)malloc(3 * sizeof(double));
	}

	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct + i) = *(edisp + 6 * n + 3 + i);
		}
		theta = vectorlength(rvct,3);
        thetaspin = spinmtx(rvct);
		thetaspin2 = matrixmatrix(thetaspin, thetaspin, 3);

		if (theta < PI / 30.0)
		{
			eta = 1.0 / 12.0 + pow(theta, 2) / 720.0 + pow(theta, 4) / 30240.0 + pow(theta, 6) / 1209600.0;
			mu = 1.0 / 360.0 + pow(theta, 2) / 7560.0 + pow(theta, 4) / 201600.0 + pow(theta, 6) / 5987520.0;
		}
		else
		{
			eta = (1.0 - 0.5 * theta / tan(0.5 * theta)) / pow(theta, 2);
			mu = (theta * theta + 4.0 * cos(theta) + theta * sin(theta) - 4.0)
				/ (4.0 * pow(theta, 4) * sin(0.5 * theta) * sin(0.5 * theta));
		}

		Ha = (double**)malloc(3 * sizeof(double*));
		for (i = 0; i < 3; i++)
		{
			*(Ha + i) = (double*)malloc(3 * sizeof(double));
		}
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(Ha + i) + j) = -0.5 * *(*(thetaspin + i) + j) + eta * *(*(thetaspin2 + i) + j);
				if (i == j)*(*(Ha + i) + j) += 1.0;
			}
		}
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(H + 6 * n + 3 + i) + 6 * n + 3 + j) = *(*(Ha + i) + j);
				if (i == j)*(*(H + 6 * n + i) + 6 * n + j) = 1.0;
			}
		}


		if(estress!=NULL && M!=NULL)
		{
			for (i = 0; i < 3; i++)
			{
				*(mvct + i) = *(estress + 6 * n + 3 + i);
			}
			dot = dotproduct(rvct,mvct,3);
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(mtheta + i) + j) = *(mvct + i) * *(rvct + j);
				}
			}
			mspin = spinmtx(mvct);
			mtheta2 = matrixmatrix(thetaspin2, mtheta, 3);
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(mtheta2 + i) + j) *= mu;
					*(*(mtheta2 + i) + j) += eta * ( *(*(mtheta + j) + i) - 2.0 * *(*(mtheta + i) + j) ) - 0.5 * *(*(mspin + i) + j);
					if (i == j)*(*(mtheta2 + i) + j) += eta * dot;
				}
			}
			Ma = matrixmatrix(mtheta2, Ha, 3);
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(M + 6 * n + 3 + i) + 6 * n + 3 + j) = *(*(Ma + i) + j);
				}
			}

			freematrix(mspin, 3);
			freematrix(mtheta2, 3);
			freematrix(Ma, 3);
		}

		freematrix(thetaspin, 3);
		freematrix(thetaspin2, 3);
		freematrix(Ha, 3);
	}

	free(rvct);
	if(estress!=NULL && M!=NULL)
	{
		free(mvct);
		freematrix(mtheta, 3);
	}
	return H;
}


double** transmatrixHPT(double* eform, double* edisp, double** T, int nnod)
{
	double** G, ** P, ** H;
	double** HP, ** HPT;

	G = spinfittermtx(eform, nnod);                     						/*SPIN-FITTER MATRIX[G].*/
	P = projectionmtx(eform, G, nnod);    										/*PROJECTION MATRIX[P].*/
	H = blockjacobimtx(edisp, NULL, NULL, nnod);								/*JACOBIAN MATRIX OF ROTATION[H].*/

	HP = matrixmatrix(H, P, 6*nnod);
	HPT = matrixmatrix(HP, T, 6*nnod);

	freematrix(G, 3);
	freematrix(P, 6*nnod);
	freematrix(H, 6*nnod);
	freematrix(HP, 6*nnod);

	return HPT;
}


double** assemgmtxCR(double* eform, double* edisp, double* estress, double* gstress, double** T, double** HPT, int nnod)
{
	int i, j, n;
	double** G, ** P, ** H;
	double** HP, ** PtHt;
	double** TtPtHt;
	double* pstress;
	double* nm, ** spinnm;
	double** Fnm, ** Fn, ** FnG, ** GtFnt;
	double** Kg, ** Kgr, ** Kgp, ** Kgm;


	Kg = (double**)malloc(6*nnod * sizeof(double*));                            /*[Kg]*/
	for (i = 0; i < 6*nnod; i++)
	{
		*(Kg + i) = (double*)malloc(6*nnod * sizeof(double));
	}

	Kgm = (double**)malloc(6*nnod * sizeof(double*));                           /*[M]&[Kgm]*/
	for (i = 0; i < 6*nnod; i++)
	{
		*(Kgm + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Kgm + i) + j) = 0.0;
		}
	}

	G = spinfittermtx(eform, nnod);                     						/*SPIN-FITTER MATRIX[G].*/
	P = projectionmtx(eform, G, nnod);    										/*PROJECTION MATRIX[P].*/
	H = blockjacobimtx(edisp, estress, Kgm, nnod);								/*JACOBIAN MATRIX OF ROTATION[H].*/

	HP = matrixmatrix(H, P, 6*nnod);                    						/*[H][P]*/
	PtHt = matrixtranspose(HP, 6*nnod);
	pstress = matrixvector(PtHt, estress, 6*nnod);     							/*projected estress {Fp}*/

	if(HPT!=NULL)
	{
		matrixmatrixII(HPT, HP, T, 6*nnod);                 						/*[H][P][T]*/
		if(gstress!=NULL)
		{
			TtPtHt = matrixtranspose(HPT, 6*nnod);
			matrixvectorII(gstress, TtPtHt, estress, 6*nnod);      						/*global estress {Fg}*/
			freematrix(TtPtHt, 6*nnod);
		}
	}


	nm = (double*)malloc(3 * sizeof(double));           						/*projected estress of each node {n}&{m}*/

	Fnm = (double**)malloc(6*nnod * sizeof(double*));                          	/*[Fnm]*/
	for (i = 0; i < 6*nnod; i++)
	{
		*(Fnm + i) = (double*)malloc(3 * sizeof(double));
	}
	Fn = (double**)malloc(6*nnod * sizeof(double*));                            /*[Fn]*/
	for (i = 0; i < 6*nnod; i++)
	{
		*(Fn + i) = (double*)malloc(3 * sizeof(double));
	}

	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(nm + i) = *(pstress + 6 * n + i);
		}
		spinnm = spinmtx(nm);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(Fnm + 6 * n + i) + j) = *(*(spinnm + i) + j);
				*(*(Fn + 6 * n + i) + j) = *(*(spinnm + i) + j);
			}
		}
		freematrix(spinnm, 3);

		for (i = 0; i < 3; i++)
		{
			*(nm + i) = *(pstress + 6 * n + 3 + i);
		}
		spinnm = spinmtx(nm);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(Fnm + 6 * n + 3 + i) + j) = *(*(spinnm + i) + j);
				*(*(Fn + 6 * n + 3 + i) + j) = 0.0;
			}
		}
		freematrix(spinnm, 3);
	}

	Kgr = matrixmatrixIII(Fnm, G, 6*nnod, 3, 6*nnod);/*[Fnm][G]*/

	FnG = matrixmatrixIII(Fn, G, 6*nnod, 3, 6*nnod);/*[Fn][G]*/
	GtFnt = matrixtranspose(FnG, 6*nnod);/*[Gt][Fnt]*/
	Kgp = matrixmatrix(GtFnt, P, 6*nnod);/*[Gt][Fnt][P]*/

	Kgm = transformationIII(Kgm, P, 6*nnod);/*[Pt][M][P]*/

	for (i = 0; i < 6*nnod; i++)
	{
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Kg + i) + j) = - *(*(Kgr + i) + j) - *(*(Kgp + i) + j) + *(*(Kgm + i) + j);
		}
	}

	Kg = transformationIII(Kg, T, 6*nnod);

	free(pstress);
	free(nm);

	freematrix(G, 3);
	freematrix(P, 6*nnod);
	freematrix(H, 6*nnod);
	freematrix(HP, 6*nnod);
	freematrix(PtHt, 6*nnod);
	freematrix(Fnm, 6*nnod);
	freematrix(Fn, 6*nnod);

	freematrix(Kgr, 6*nnod);
	freematrix(FnG, 6*nnod);
	freematrix(GtFnt, 6*nnod);
	freematrix(Kgp, 6*nnod);
	freematrix(Kgm, 6*nnod);
	return Kg;
}


double** assemtmtxCR(double** Ke, double* eform, double* edisp, double* estress, double* gstress, double** T, double** HPT, int nnod)
{
	int i, j;
	double** Kt;

	Kt = assemgmtxCR(eform, edisp, estress, gstress, T, HPT, nnod);/*[Kg]=[Kgr]+[Kgp]+[Kgm]*/
	Ke = transformationIII(Ke, HPT, 6*nnod);/*[Ke]=[Tt][Pt][Ht][K][H][P][T]*/

	for (i = 0; i < 6*nnod; i++)
	{
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Kt + i) + j) += *(*(Ke + i) + j);/*[Kt]=[Ke]+[Kg]*/
		}
	}
	return Kt;
}


void symmetricmtx(double** estiff, int msize)
{
	int i, j;
	for (i = 0; i < msize; i++)
	{
		for (j = 0; j < i; j++)
		{
			*(*(estiff + i) + j) = 0.5 * (*(*(estiff + i) + j) + *(*(estiff + j) + i));
			*(*(estiff + j) + i) = *(*(estiff + i) + j);
		}
	}
	return;
}


void updaterotation(double* ddisp, double* gvct, int nnode)
/*FORMATION UPDATE IF ROTATION IS FINITE.*/
{
	int i, j;
	long int loff;
	double* rvctR, * rvctL, * rvct;
	double** rmtxR, ** rmtxL, ** rmtx;
	rvctR = (double*)malloc(3 * sizeof(double));
	rvctL = (double*)malloc(3 * sizeof(double));
	for (i = 0; i < nnode; i++)
	{
		for (j = 0; j < 3; j++)
		{
			loff = 6 * i + j;
			*(ddisp + loff) += *(gvct + loff);
		}
		for (j = 0; j < 3; j++)
		{
			loff = 6 * i + 3 + j;
			*(rvctR + j) = *(ddisp + loff);
			*(rvctL + j) = *(gvct + loff);
		}
		rmtxR = rotationmtx(rvctR);
		rmtxL = rotationmtx(rvctL);
		rmtx = matrixmatrix(rmtxL, rmtxR, 3);
		rvct = rotationvct(rmtx);
		for (j = 0; j < 3; j++)
		{
			loff = 6 * i + 3 + j;
			*(ddisp + loff) = *(rvct + j);
		}
		freematrix(rmtxR, 3);
		freematrix(rmtxL, 3);
		freematrix(rmtx, 3);
		free(rvct);
	}
	free(rvctR);
	free(rvctL);
	return;
}/*updaterotation*/



double* quaternionvct(double* rvct)
{
	int i;
	double* qvct;
	double theta;

	qvct = (double*)malloc(4 * sizeof(double));
	theta = sqrt(*(rvct + 0) * *(rvct + 0) + *(rvct + 1) * *(rvct + 1) + *(rvct + 2) * *(rvct + 2));
	if(theta!=0)
	{
		for (i = 0; i < 3; i++)
		{
			*(qvct + i) = sin(0.5*theta)**(rvct + i)/theta;
		}
	}
	else
	{
        for (i = 0; i < 3; i++)
		{
			*(qvct + i) = 0.0;
		}
    }
	*(qvct + 3) = cos(0.5*theta);
	return qvct;
}


double** updatedrccos(double** drccosinit, double* gforminit, double* gform)
{
	int i;
	double** drccos;
	double* rvct1, * rvct2;
	double** rmtxinit, ** rmtxinitt, **rmtx, ** drccosrmtx;

	rvct1 = (double*)malloc(3 * sizeof(double));
	rvct2 = (double*)malloc(3 * sizeof(double));

	for (i = 0; i < 3; i++)
	{
		*(rvct1 + i) = *(gforminit + 3 + i);
		*(rvct2 + i) = *(gforminit + 9 + i);
	}
	rmtxinit = interpolatermtx(rvct1, rvct2, 0.5);
	rmtxinitt = matrixtranspose(rmtxinit, 3);

	for (i = 0; i < 3; i++)
	{
		*(rvct1 + i) = *(gform + 3 + i);
		*(rvct2 + i) = *(gform + 9 + i);
	}
	rmtx = interpolatermtx(rvct1, rvct2, 0.5);

	drccosrmtx = matrixmatrix(rmtx, rmtxinitt, 3);
	drccos = matrixmatrix(drccosrmtx, drccosinit, 3);

	free(rvct1);
	free(rvct2);

	freematrix(rmtxinit, 3);
	freematrix(rmtxinitt, 3);
	freematrix(rmtx, 3);
	freematrix(drccosrmtx, 3);

	return drccos;
}


double** interpolatermtx(double* rvct1, double* rvct2, double alpha)
{
	int i;
	double* rvct;
	double** rmtx1, ** rmtx2, **trmtx1, ** rmtx;
	double** alpharmtx, ** midrmtx;


	/*
	qvct1 = quaternionvct(rvct1);
	qvct2 = quaternionvct(rvct2);
	dot = dotproduct(qvct1,qvct2,4);
	theta = acos(dot);

	if(theta!= 0)
	{
		for (i = 0; i < 4; i++)
		{
			*(midqvct + i) = (*(qvct1 + i)+*(qvct2 + i))/(2*cos(0.5*theta));
		}
	}
	*/



	rmtx1 = rotationmtx(rvct1);
	rmtx2 = rotationmtx(rvct2);
	trmtx1 = matrixtranspose(rmtx1, 3);
	rmtx = matrixmatrix(rmtx2, trmtx1, 3);
	rvct = rotationvct(rmtx);

	for (i = 0; i < 3; i++)
	{
		*(rvct+i)*=alpha;
	}

	alpharmtx = rotationmtx(rvct);
	midrmtx = matrixmatrix(alpharmtx, rmtx1, 3);


	freematrix(rmtx1, 3);
	freematrix(rmtx2, 3);
	freematrix(trmtx1, 3);
	freematrix(rmtx, 3);
	free(rvct);
	freematrix(alpharmtx, 3);

	return midrmtx;
}



double* extractlocalcoord(double* gform, double** drccos, double nnod)
/*EXTRACT LOCAL ELEMENT DEFORMATION FROM GLOBAL.*/
/*UPDATE PSUEDO-ROTATION VECTOR*/
{
	long int i, n;
	double* d, * r, * c, * td, * tr, * eform;
	double** trmtx, ** rmtx;

	eform = (double*)malloc(6 * nnod * sizeof(double));

	c = (double*)malloc(3 * sizeof(double));
	d = (double*)malloc(3 * sizeof(double));
	r = (double*)malloc(3 * sizeof(double));

	for (i = 0; i < 3; i++)
	{
		*(c + i) = 0.0;
		for (n = 0; n < nnod; n++)
		{
			*(c + i) += *(gform + 6 * n + i) / nnod;
		}
	}
	/*CENTER*/
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(d + i) = *(gform + 6 * n + i) - *(c + i);
			*(r + i) = *(gform + 6 * n + 3 + i);
		}

		td = matrixvector(drccos, d, 3);
		/*EACH NODE FROM CENTER*/

		rmtx = rotationmtx(r);

		/*rmtx:Ra*/
		trmtx = matrixmatrix(drccos, rmtx, 3);/*TRIAD DIRECTION MATRIX(3 VECTOR) IN LOCAL*/
		tr = rotationvct(trmtx);

		for (i = 0; i < 3; i++)
		{
			*(eform + 6 * n + i) = *(td + i);
			*(eform + 6 * n + 3 + i) = *(tr + i);
		}
		free(td);
		free(tr);
		freematrix(rmtx, 3);
		freematrix(trmtx, 3);
	}
	free(c);
	free(d);
	free(r);

	return eform;
}/*extractlocalcoord*/



double*  extractdeformation(double* eforminit, double* eform, int nnod)
/*EXTRACT LOCAL ELEMENT DEFORMATION FROM GLOBAL.*/
/*UPDATE PSUEDO-ROTATION VECTOR*/
{
	int n, i;
	double* edisp;
	double* r, * rinit, * rvct;
	double** rmtx, ** rh, ** rt, ** rtt;

	edisp = (double*)malloc(6*nnod * sizeof(double));
	r     = (double*)malloc(3 * sizeof(double));
	rinit = (double*)malloc(3 * sizeof(double));

	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(r + i)     = *(eform     + 6 * n + 3 + i);
			*(rinit + i) = *(eforminit + 6 * n + 3 + i);
		}

		rh = rotationmtx(r);
		rt = rotationmtx(rinit);
		rtt = matrixtranspose(rt, 3);
		rmtx = matrixmatrix(rh, rtt, 3);
		rvct = rotationvct(rmtx);

		for (i = 0; i < 3; i++)
		{
			*(edisp + 6 * n + i)     = *(eform + 6 * n + i) - *(eforminit + 6 * n + i);
			*(edisp + 6 * n + 3 + i) = *(rvct + i);
		}

		freematrix(rh, 3);
		freematrix(rt, 3);
		freematrix(rtt, 3);
		freematrix(rmtx, 3);
		free(rvct);
	}

	free(r);
	free(rinit);
	return edisp;
}



/*MATERIAL & SPATIAL FORM VARIABLES.*/


double* pullback(double* ddisp, double* gvct_s, int nnode)
{
	int i,n;
	double* rvct, * vct_s, * vct_m, * gvct_m;
	double** rmtx, ** trmtx;

	gvct_m = (double*)malloc(6*nnode * sizeof(double));
	vct_s = (double*)malloc(3 * sizeof(double));
	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnode; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct + i)  = *(ddisp + 6*n+3+i);
			*(vct_s + i) = *(gvct_s + 6*n+3+i);
		}
		rmtx = rotationmtx(rvct);
		trmtx = matrixtranspose(rmtx, 3);
		vct_m = matrixvector(trmtx, vct_s, 3);
		for (i = 0; i < 3; i++)
		{
			*(gvct_m + 6*n+i)   = *(gvct_s + 6*n+i);
			*(gvct_m + 6*n+3+i) = *(vct_m + i);
		}
		freematrix(rmtx,3);
		freematrix(trmtx,3);
		free(vct_m);
	}
	free(rvct);
	free(vct_s);
	return gvct_m;
}

double** pullbackmtx(double* gform, int nnod)
{
	int i,j,n;
	double* rvct;
	double** rmtx,** trmtx, **Rt;

	Rt = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(Rt + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Rt + i) + j) = 0.0;
		}
	}

	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct  + i) = *(gform + 6*n+3+i);
		}
		rmtx = rotationmtx(rvct);
		trmtx = matrixtranspose(rmtx, 3);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(Rt + 6 * n + 3 + i) + 6 * n + 3 + j) = *(*(trmtx + i) + j);
				if (i == j)*(*(Rt + 6 * n + i) + 6 * n + j) = 1.0;
			}
		}
		freematrix(rmtx,3);
		freematrix(trmtx,3);
	}
	free(rvct);
	return Rt;
}

double* pushforward(double* ddisp, double* gvct_m, int nnode)
{
	int i,n;
	double* rvct, * vct_s, * vct_m, * gvct_s;
	double** rmtx;

	gvct_s = (double*)malloc(6*nnode * sizeof(double));
	vct_m = (double*)malloc(3 * sizeof(double));
	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnode; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct  + i) = *(ddisp + 6*n+3+i);
			*(vct_m + i) = *(gvct_m + 6*n+3+i);
		}
		rmtx = rotationmtx(rvct);
		vct_s = matrixvector(rmtx, vct_m, 3);
		for (i = 0; i < 3; i++)
		{
			*(gvct_s + 6*n+i) = *(gvct_m + 6*n+i);
			*(gvct_s + 6*n+i+3) = *(vct_s + i);
		}
		free(vct_s);
		freematrix(rmtx,3);
	}
	free(rvct);
	free(vct_m);
	return gvct_s;
}


double** pushforwardmtx(double* gform, int nnod)
{
	int i,j,n;
	double* rvct;
	double** rmtx, **R;

	R = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(R + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(R + i) + j) = 0.0;
		}
	}
	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct  + i) = *(gform + 6*n+3+i);
		}
		rmtx = rotationmtx(rvct);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(R + 6 * n + 3 + i) + 6 * n + 3 + j) = *(*(rmtx + i) + j);
				if (i == j)*(*(R + 6 * n + i) + 6 * n + j) = 1.0;
			}
		}
		freematrix(rmtx,3);
	}
	free(rvct);

	return R;
}

/*MID-POINT VARIABLES.*/
double* midpointvct(double* vct,double* lastvct,double alpha,int size)
{
	int i;
	double* midvct;

	midvct = (double*)malloc(size * sizeof(double));
	for (i=0;i<size;i++)
	{
		*(midvct+i)=(1.0-alpha)**(vct+i)+alpha**(lastvct+i);
	}
	return midvct;
}

double** midpointmtx(double** mtx,double** lastmtx,double alpha,int size)
{
	int i,j;
	double** midmtx;

	midmtx= (double**)malloc(size * sizeof(double*));
	for (i=0;i<size;i++)
	{
		*(midmtx+i) = (double*)malloc(size * sizeof(double));
		for(j=0;j<size;j++)
		{
			*(*(midmtx+i)+j)=(1.0-alpha)**(*(mtx+i)+j) + alpha**(*(lastmtx+i)+j);
		}
	}
	return midmtx;
}


