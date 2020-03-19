
#include <iostream>

#include "problemmodel.h"
#include "mathfunction.h"

void PoissonProblem::Source(double x,double y,double *value)
{
	*value=2*PI*PI*sin(PI*x)*sin(PI*y);
//	value[0]=-2*(y*y-y+x*x-x);
}

void PoissonProblem::DirichletBoundary(double x,double y,double *value)
{
	value[0]=0;
}
void PoissonProblem::Solution(double x,double y,double *value)
{
	value[0]=sin(PI*x)*sin(PI*y);
	value[1]=PI*cos(PI*x)*sin(PI*y);
	value[2]=PI*sin(PI*x)*cos(PI*y);
/*	value[0]=x*y*(x-1)*(y-1);
	value[1]=(2*x-1)*y*(y-1);
	value[2]=x*(x-1)*(2*y-1);
*/
}

void PoissonModel::GetLocRHS(int ElemID,double *LocF)
{
	VE.GetRHS(ElemID,pde.Source,LocF);
}

void PoissonModel::Solve()
{   
    GetSystem();
    //deal with boundary condition
	int i,k,BdofNum=ms.B_Node_Num*dof.Num_PerNode+ms.B_Edge_Num*dof.Num_PerEdge,IdofNum=dof.Dof_Num-BdofNum;
	int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dof.Dof_Num];
	double *BdofVal=new double[BdofNum];
	for(i=0;i<BdofNum;i++) BdofVal[i]=0;
	VE.GetBdof_BdofVal(pde.DirichletBoundary,Bdof,BdofVal);
 
    long N=StiffMatrix.nonZeros();
    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(N);
    
    for(i=0;i<BdofNum;i++)
        uh[Bdof[i]]=BdofVal[i];
 
    RHS-=StiffMatrix*uh.pVector;
  
    VectorXd b(IdofNum),x(IdofNum);
    int col=0,row;
    for(i=0;i<dof.Dof_Num;i++) //get the relation between dof and inner dof
    {
        IdofID[i]=-1;
        if(dof.DType[i]==0)
        {
            b[col]=RHS[i];
            Idof[col]=i;
            IdofID[i]=col;
            col++;
        }
    }
    
    for(col=0;col<IdofNum;col++)
    {
            i=Idof[col];
            for(SparseMatrix<double>::InnerIterator it(StiffMatrix,i);it;++it)
            {
                k=(int) it.row();
                if(dof.DType[k]==0)
                {
                    row=IdofID[k];
                    tripletList.push_back(Eigen::Triplet<double>(row,col,it.value()));
                }
            }
    }
    
    SparseMatrix<double> A(IdofNum,IdofNum);
    A.setFromTriplets(tripletList.begin(),tripletList.end());

    //direct method
    Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    
    // iterative method
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    //    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    x=solver.solve(b);
    for(i=0;i<IdofNum;i++) uh[Idof[i]]=x[i];
    
    delete []Bdof; delete []BdofVal; delete []Idof; delete []IdofID;
}

void PoissonModel::GetSystem()
{
	int i,j,k,EdofNum,N=0;
    for(i=0;i<ms.Element_Num;i++)
    {
        j=VE.dof.Total_Num_PerElement[i];  N+=j*j;
    }
    std::vector< Eigen::Triplet<double> > tripletList; //for constructing sparsmatrix
    tripletList.reserve(N);
    
	for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
	{
		EdofNum=dof.Total_Num_PerElement[i];
		double **LocA=new double*[EdofNum];
		double *LocF=new double[EdofNum];
		for(j=0;j<EdofNum;j++)
			LocA[j]=new double[EdofNum];
		VE.GetA_H1(i,LocA);
		GetLocRHS(i,LocF);
		for(j=0;j<EdofNum;j++)//组装
		{
			RHS[ dof.TotalD[i][j] ]+=LocF[j];//组装右端项
			for(k=0;k<EdofNum;k++)
                tripletList.push_back(Eigen::Triplet<double>(dof.TotalD[i][j],dof.TotalD[i][k],LocA[j][k]));
		}
		for(j=0;j<EdofNum;j++)
			delete [] LocA[j];
		delete []LocA; delete []LocF;			
	}
    StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
}

double PoissonModel::EnergyError()
{
    VectorXd X=uI.pVector-uh.pVector;
	double value=X.transpose()*StiffMatrix*X;
	return sqrt(value);
}

double PoissonModel::MaxError()
{
	double value=0;
	int i;
	for(i=0;i<dof.Dof_Num;i++)
		if(value<fabs(uI[i]-uh[i]))
			value=fabs(uI[i]-uh[i]);
	return value;
}

/*another method to deal boundary conditor
 
 void PoissonModel::GetSystem()
 {
 int i,j,k,EdofNum,N=dof.Dof_NumofAllElements;
 
 // begin to prepare for dealing with boundary condition
 int BdofNum=ms.B_Node_Num*dof.Num_PerNode+ms.B_Edge_Num*dof.Num_PerEdge,IdofNum=dof.Dof_Num-BdofNum;
 int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dof.Dof_Num];
 double *BdofVal=new double[BdofNum];//Idof: inner dof to dof ID, IdofID: dof to inner dof ID
 for(i=0;i<BdofNum;i++) BdofVal[i]=0;
 int col=0,row;
 for(i=0;i<dof.Dof_Num;i++) //get the relation between dof and inner dof
 {
 IdofID[i]=-1;
 if(dof.DType[i]==0)
 {
 Idof[col]=i;  IdofID[i]=col; col++;
 }
 }
 // end of preparation
 
 std::vector< Eigen::Triplet<double> > tripletList; //for constructing sparsmatrix
 tripletList.reserve(N);
 
 std::vector< Eigen::Triplet<double> > tripletList2; //for constructing sparsmatrix wrt. inner dof
 tripletList2.reserve(N);
 
 for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
 {
 EdofNum=dof.Total_Num_PerElement[i];
 double **LocA=new double*[EdofNum];
 double *LocF=new double[EdofNum];
 for(j=0;j<EdofNum;j++)
 LocA[j]=new double[EdofNum];
 VE.GetA_H1(i,LocA);
 GetLocRHS(i,LocF);
 for(j=0;j<EdofNum;j++)//组装
 {
 row= dof.TotalD[i][j];
 RHS[ row ]+=LocF[j];//组装右端项
 for(k=0;k<EdofNum;k++)
 {
 col=dof.TotalD[i][k] ;
 tripletList.push_back(Eigen::Triplet<double>(row,col,LocA[j][k]));
 if(dof.DType[row]==0&&dof.DType[col]==0)
 tripletList2.push_back(Eigen::Triplet<double>(IdofID[row],IdofID[col],LocA[j][k]));
 }
 }
 
 for(j=0;j<EdofNum;j++)
 delete [] LocA[j];
 delete []LocA; delete []LocF;
 }
 StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
 
 //deal with boundary condition
 VE.GetBdof_BdofVal(pde.DirichletBoundary,Bdof,BdofVal);
 for(i=0;i<BdofNum;i++)
 uh[Bdof[i]]=BdofVal[i];
 RHS-=StiffMatrix*uh.pVector;
 VectorXd b(IdofNum),x(IdofNum);
 b.setZero();
 for(i=0;i<dof.Dof_Num;i++)
 if(dof.DType[i]==0)
 b[ IdofID[i] ]=RHS[i];
 
 SparseMatrix<double> A(IdofNum,IdofNum);
 A.setFromTriplets(tripletList2.begin(),tripletList2.end());
 
 //direct method
 Eigen::SimplicialLLT<SparseMatrix<double> > solver;
 //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
 
 // iterative method
 //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
 //    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
 //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
 solver.compute(A);
 x=solver.solve(b);
 for(i=0;i<IdofNum;i++) uh[Idof[i]]=x[i];
 
 delete []Bdof; delete []BdofVal; delete []Idof; delete []IdofID;
 
 }
 
 */
