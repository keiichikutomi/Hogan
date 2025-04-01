
double* assemeinternal(struct owire* elem)
{
  int ii,jj,i,j;
  int nnod;
  double* einternal,*gpstress,*gpinternal;
  double** Bt;

  nnod = elem->nnod;


  einternal = (double*)malloc(6 * nnod * sizeof(double));
  for(i = 0; i < nnod; i++)
  {
	for(j = 0; j < 6; j++)
	{
	   *(einternal+6*i+j) = 0.0/*shell->stress[i][j]*/;
	}
  }



  for(i = 0; i < nnod; i++)
  {
	for(j = 0; j < 6; j++)
	{
	  elem->stress[i][j] = *(einternal+6*i+j);
	}
  }
  return einternal;
}

