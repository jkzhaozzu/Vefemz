#include "stdlib.h"
#include "time.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers> 
#include <Eigen/SparseQR> 
#include "mathfunction.h"
#include "vemfunction.h"
/*
void toEigenSparseSolve(Matrix &StiffMatrix,Vector &RHS,VEMFunction &uh)
{
	
	int N=StiffMatrix.GetRealSize();
	std::vector< Eigen::Triplet<double> > tripletList;
	tripletList.reserve(N);
	M_Element *p;
	for(int i=0;i<StiffMatrix.size;i++)
    {
        p=StiffMatrix.Header[i].pNext;
        while(p!=NULL)
        {
			tripletList.push_back(Eigen::Triplet<double>(i,p->pos,p->value));
            p=p->pNext;
        }
    }
	
	Eigen::VectorXd x(uh.dof.Dof_Num), b(uh.dof.Dof_Num);
	Eigen::SparseMatrix<double> A(uh.dof.Dof_Num,uh.dof.Dof_Num);
	A.setFromTriplets(tripletList.begin(),tripletList.end());
	for(int i=0;i<uh.dof.Dof_Num;i++)
		b[i]=RHS[i];
//direct method
//	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

// iterative method
//	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
//	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
//	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	solver.compute(A);
	x = solver.solve(b);

	for(int i=0;i<uh.dof.Dof_Num;i++)
		uh[i]=x[i];
}

void toEigenSparseSolve(Matrix &StiffMatrix,Vector &RHS,VEMFunction &uh,VEMFunction &ph)
{
	
	int N=StiffMatrix.GetRealSize();
	std::vector< Eigen::Triplet<double> > tripletList;
	tripletList.reserve(N);
	M_Element *p;
	for(int i=0;i<StiffMatrix.size;i++)
    {
        p=StiffMatrix.Header[i].pNext;
        while(p!=NULL)
        {
			tripletList.push_back(Eigen::Triplet<double>(i,p->pos,p->value));
            p=p->pNext;
        }
    }
	
	Eigen::VectorXd x(StiffMatrix.size), b(StiffMatrix.size);
	Eigen::SparseMatrix<double> A(StiffMatrix.size,StiffMatrix.size);
	A.setFromTriplets(tripletList.begin(),tripletList.end());
	for(int i=0;i<StiffMatrix.size;i++)
		b[i]=RHS[i];
//direct method
//	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
//	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

// iterative method
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
//	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
//	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	solver.compute(A);
	x = solver.solve(b);

	for(int i=0;i<uh.dof.Dof_Num;i++)
		uh[i]=x[i];
	for(int i=uh.dof.Dof_Num;i<StiffMatrix.size;i++)
		ph[i-uh.dof.Dof_Num]=x[i];

}

void toEigenSparseSolveL02(Matrix &StiffMatrix,Vector &RHS,VEMFunction &uh,VEMFunction &ph)
{
	
	int N=StiffMatrix.GetRealSize();
	std::vector< Eigen::Triplet<double> > tripletList;
	tripletList.reserve(N+ph.ms.Element_Num);
	M_Element *p;
	for(int i=0;i<StiffMatrix.size;i++)
    {
        p=StiffMatrix.Header[i].pNext;
        while(p!=NULL)
        {
			tripletList.push_back(Eigen::Triplet<double>(i,p->pos,p->value));
            p=p->pNext;
        }
    }
	for(int i=0;i<ph.ms.Element_Num;i++)
		tripletList.push_back(Eigen::Triplet<double>(StiffMatrix.size,uh.dof.Dof_Num+i*ph.dof.Num_PerElement,ph.ms.ElementMeasure[i]));
	
	Eigen::VectorXd x(StiffMatrix.size), b(StiffMatrix.size+1);
	Eigen::SparseMatrix<double> A(StiffMatrix.size+1,StiffMatrix.size);
	A.setFromTriplets(tripletList.begin(),tripletList.end());
	for(int i=0;i<StiffMatrix.size;i++)
		b[i]=RHS[i];
	b[StiffMatrix.size]=0;
//direct method
//	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
//	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

// iterative method
//	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
//	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	solver.compute(A);
	x = solver.solve(b);

	for(int i=0;i<uh.dof.Dof_Num;i++)
		uh[i]=x[i];
	for(int i=uh.dof.Dof_Num;i<StiffMatrix.size;i++)
		ph[i-uh.dof.Dof_Num]=x[i];

}
 */
