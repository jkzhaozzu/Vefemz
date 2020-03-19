//
//  FEpoissonmodel.cpp
//  Vefemz
//
//  Created by 张蓓 on 2019/11/10.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//
#include <iostream>
#include "FEmodel.h"

void FEPoissonModel::GetLocRHS(int ElemID,double *LocF)
{
    FE.GetRHS(ElemID,pde.Source,LocF);
}

void FEPoissonModel::GetSystem()
{
    int i,j,k,EdofNum=FE.dof.Total_Num_PerElement[0],N=ms.Element_Num*EdofNum*EdofNum;
    std::vector< Eigen::Triplet<double> > tripletList; //for constructing sparsmatrix
    tripletList.reserve(N);
    
    for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
    {
        double **LocA=new double*[EdofNum];
        double *LocF=new double[EdofNum];
        for(j=0;j<EdofNum;j++)
            LocA[j]=new double[EdofNum];
        FE.GetA_H1(i,LocA);
        GetLocRHS(i,LocF);
        for(j=0;j<EdofNum;j++)//组装
        {
            RHS[ dof.TotalD[i][j] ]+=LocF[j];//组装右端项
            for(k=0;k<EdofNum;k++)
                tripletList.push_back(Eigen::Triplet<double>(dof.TotalD[i][j],dof.TotalD[i][k],LocA[j][k]));
 //           for(k=0;k<EdofNum;k++)std::cout<<LocA[j][k]<<std::endl;
        }
        for(j=0;j<EdofNum;j++)
            delete [] LocA[j];
        delete []LocA; delete []LocF;
    }
    StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
}

void FEPoissonModel::Solve()
{
    GetSystem();
    //deal with boundary condition
    int i,k,BdofNum=ms.B_Node_Num*dof.Num_PerNode+ms.B_Edge_Num*dof.Num_PerEdge,IdofNum=dof.Dof_Num-BdofNum;
    int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dof.Dof_Num];
    double *BdofVal=new double[BdofNum];
    for(i=0;i<BdofNum;i++) BdofVal[i]=0;
    FE.GetBdof_BdofVal(pde.DirichletBoundary,Bdof,BdofVal);
    
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

double FEPoissonModel::EnergyError()
{
    VectorXd X=uI.pVector-uh.pVector;
    double value=X.transpose()*StiffMatrix*X;
    return sqrt(value);
}

double FEPoissonModel::MaxError()
{
    double value=0;
    int i;
    for(i=0;i<dof.Dof_Num;i++)
        if(value<fabs(uI[i]-uh[i]))
            value=fabs(uI[i]-uh[i]);
    return value;
}


