#include <iostream>
#include <Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>

#include "problemmodel.h"
#include "mathfunction.h"
#include "sparsesolve.h"


void ParabolicProblem::Source(double x,double y,double t,double *value)
{
    *value=exp(t)*sin(PI*x)*sin(PI*y)+2*PI*PI*exp(t)*sin(PI*x)*sin(PI*y);
    //    *value=10*exp(10*t)*sin(PI*x)*sin(PI*y)+2*PI*PI*exp(10*t)*sin(PI*x)*sin(PI*y);
}

void ParabolicProblem::DirichletBoundary(double x,double y,double t,double *value)
{
    value[0]=0;
}
void ParabolicProblem::InitialData(double x,double y,double *value)
{
    value[0]=sin(PI*x)*sin(PI*y);
    value[1]=PI*cos(PI*x)*sin(PI*y);
    value[2]=PI*sin(PI*x)*cos(PI*y);
    
    /*    value[0]=sin(PI*x)*sin(PI*y);
     value[1]=PI*cos(PI*x)*sin(PI*y);
     value[2]=PI*sin(PI*x)*cos(PI*y);
     */
}
void ParabolicProblem::Solution(double x,double y,double t,double *value)
{
    value[0]=exp(t)*sin(PI*x)*sin(PI*y);
    value[1]=PI*exp(t)*cos(PI*x)*sin(PI*y);
    value[2]=PI*exp(t)*sin(PI*x)*cos(PI*y);
    
    /*    value[0]=exp(10*t)*sin(PI*x)*sin(PI*y);
     value[1]=PI*exp(10*t)*cos(PI*x)*sin(PI*y);
     value[2]=PI*exp(10*t)*sin(PI*x)*cos(PI*y);
     */
}



void ParabolicModel::GetLocRHS(int ElemID,double t,double *LocF)
{
    
    VE.GetRHS(ElemID,pde.Source,t,LocF);
}

void ParabolicModel::GetMatrix()
{
    int i,j,k,EdofNum,N=0;
    for(i=0;i<ms.Element_Num;i++)
    {
        j=VE.dof.Total_Num_PerElement[i];  N+=j*j;
    }
    std::vector< Eigen::Triplet<double> > tripletList; //for constructing StiffMatrix
    tripletList.reserve(N);
    std::vector< Eigen::Triplet<double> > tripletListH1; //for constructing MatrixH1
    tripletListH1.reserve(N);
    std::vector< Eigen::Triplet<double> > tripletListL2; //for constructing MatrixL2
    tripletListL2.reserve(N);
    
    for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
    {
        EdofNum=dof.Total_Num_PerElement[i];
        double **LocAH1=new double*[EdofNum]; double **LocAL2=new double*[EdofNum];
        for(j=0;j<EdofNum;j++)
        {
            LocAH1[j]=new double[EdofNum];LocAL2[j]=new double[EdofNum];
        }
        VE.GetA_H1(i,LocAH1);    VE.GetA_L2(i,LocAL2);
        
        for(j=0;j<EdofNum;j++)//组装
            for(k=0;k<EdofNum;k++)
            {
                tripletList.push_back(Eigen::Triplet<double>(dof.TotalD[i][j],dof.TotalD[i][k],LocAH1[j][k]+LocAL2[j][k]/dt));
                tripletListH1.push_back(Eigen::Triplet<double>(dof.TotalD[i][j],dof.TotalD[i][k],LocAH1[j][k]));
                tripletListL2.push_back(Eigen::Triplet<double>(dof.TotalD[i][j],dof.TotalD[i][k],LocAL2[j][k]));
            }
        for(j=0;j<EdofNum;j++)
        {
            delete [] LocAH1[j]; delete [] LocAL2[j];
        }
        delete []LocAH1;     delete []LocAL2;
    }
    StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
    MatrixH1.setFromTriplets(tripletListH1.begin(),tripletListH1.end());
    MatrixL2.setFromTriplets(tripletListL2.begin(),tripletListL2.end());
}

void ParabolicModel::GetLV(double t)
{
    int i,j,EdofNum;
    RHS.setZero();
    for(i=0;i<ms.Element_Num;i++)
    {
        EdofNum=dof.Total_Num_PerElement[i];
        double *LocF=new double[EdofNum];
        GetLocRHS(i,t,LocF);
        for(j=0;j<EdofNum;j++)
            RHS[ dof.TotalD[i][j] ]+=LocF[j];
        delete []LocF;
    }
}

void ParabolicModel::Solve()
{
    GetMatrix();
    //deal with boundary condition
    int i,mi,k,BdofNum=ms.B_Node_Num*dof.Num_PerNode+ms.B_Edge_Num*dof.Num_PerEdge,IdofNum=dof.Dof_Num-BdofNum;
    int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dof.Dof_Num];
    double *BdofVal=new double[BdofNum];
     
    long N=StiffMatrix.nonZeros();
    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(N);

    int col=0,row;
    for(i=0;i<dof.Dof_Num;i++) //get the relation between dof and inner dof
    {
        IdofID[i]=-1;
        if(dof.DType[i]==0)
        {
            Idof[col]=i;  IdofID[i]=col; col++;
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
    VectorXd b(IdofNum),x(IdofNum);
    //direct method
    Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
     
    // iterative method
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    //    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
     
    for(mi=0;mi<dtn;mi++)//iteration solving starts
    {
        tn0=mi*dt; tn=(mi+1)*dt;
        uh.pVector.setZero();
        GetLV(tn);
        RHS+=(MatrixL2*uh0.pVector/dt);
        for(i=0;i<BdofNum;i++)
            BdofVal[i]=0;
        VE.GetBdof_BdofVal(pde.DirichletBoundary,tn,Bdof,BdofVal);
        for(i=0;i<BdofNum;i++)
            uh[Bdof[i]]=BdofVal[i];
     
        RHS-=StiffMatrix*uh.pVector;
       for(i=0;i<dof.Dof_Num;i++) //get the relation between dof and inner dof
            if(dof.DType[i]==0)
                b[ IdofID[i] ]=RHS[i];
        solver.compute(A);
        x=solver.solve(b);
        for(i=0;i<IdofNum;i++) uh[Idof[i]]=x[i];
        uh0.pVector=uh.pVector;
    }
}

double ParabolicModel::ErrorH1()
{
    VectorXd X=uI.pVector-uh.pVector;
    double value=X.transpose()*MatrixH1*X;
    return sqrt(value);
}

double ParabolicModel::ErrorL2()
{
    VectorXd X=uI.pVector-uh.pVector;
    double value=X.transpose()*MatrixL2*X;
    return sqrt(value);
}


