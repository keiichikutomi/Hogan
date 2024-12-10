

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
//#include <Eigen/Dense>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>

#include <iostream>
#include <vector>


using SparseMatrix = Eigen::SparseMatrix<double,Eigen::ColMajor,int>;
using Vector = Eigen::VectorXd;
using Triplet = Eigen::Triplet<double>;



void assemtriplet(std::vector<Triplet>& triplets,
				  double** estiff,
				  struct owire *elem,
				  long int *constraintmain,
				  struct oconf *confs)
{
  long int i, j, ii, jj, n1, n2;
  double edata;
  int nnod=elem->nnod;
  for(n1=0;n1<nnod;n1++)
  {
	for(i=0;i<6;i++)
	{
	  ii=6*(elem->node[n1]->loff)+i;
	  if(constraintmain!=NULL)ii=*(constraintmain+ii);

	  for(n2=0;n2<nnod;n2++)
	  {
		for(j=0;j<6;j++)
		{
		  jj=6*(elem->node[n2]->loff)+j;
		  if(constraintmain!=NULL)jj=*(constraintmain+jj);

		  //if(ii>=jj)
		  //{
			if((confs+ii)->iconf!=1 && (confs+jj)->iconf!=1)
			{
				edata=*(*(estiff+6*n1+i)+6*n2+j);
				if(edata!=0.0)
				{
					triplets.emplace_back(ii, jj, edata);
				}
			}
			else if(ii==jj)
			{
				triplets.emplace_back(ii, jj, 1.0);
			}
		  //}
		}
	  }
	}
  }

}

void assemtripletshell(std::vector<Triplet>& triplets,
					   double** estiff,
					   struct oshell *shell,
					   long int *constraintmain,
					   struct oconf *confs)
{
  long int i, j, ii, jj, n1, n2;
  double edata;
  int nnod=shell->nnod;
  for(n1=0;n1<nnod;n1++)
  {
	for(i=0;i<6;i++)
	{
	  ii=6*(shell->node[n1]->loff)+i;
	  if(constraintmain!=NULL)ii=*(constraintmain+ii);

	  for(n2=0;n2<nnod;n2++)
	  {
		for(j=0;j<6;j++)
		{
		  jj=6*(shell->node[n2]->loff)+j;
		  if(constraintmain!=NULL)jj=*(constraintmain+jj);

		  //if(ii>=jj)
		  //{
			if((confs+ii)->iconf!=1 && (confs+jj)->iconf!=1)
			{
				edata=*(*(estiff+6*n1+i)+6*n2+j);
				if(edata!=0.0)
				{
					triplets.emplace_back(ii, jj, edata);
				}
			}
			else if(ii==jj)
			{
				triplets.emplace_back(ii, jj, 1.0);
			}
		  //}
		}
	  }
	}
  }

}

/*
void boundary(SparseMatrix& gmtx, struct oconf* confs)
{
	for(int i=0;i<gmtx.outerSize();++i)
	{
		for(SparseMatrix::InnerIterator it(gmtx, i); it; ++it)
		{
			int id_col = it.col(); // —ñ”Ô†
			int id_row = it.row(); // s”Ô†

			if(i == it.col()){
				it.valueRef() = 1.0;
			}
			else
			{
				it.valueRef() = 0.0;      //value()‚Æ‚Ìˆá‚¢
			}
		}
	}
}
*/







/*


SparseMatrix& K_global,
Vector& F,
const std::vector<Triplet>& triplets


gmtx.coeffRef(0, 0) = 1.0;
gmtx.prune([](int i, int j, double) { return i != j; });


*/





