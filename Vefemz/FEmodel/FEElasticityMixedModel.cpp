//
//  FEElasticityMixedModel.cpp
//  Vefemz
//
//  Created by 张蓓 on 2020/2/18.
//  Copyright © 2020年 Jikun Zhao. All rights reserved.
//

#include <iostream>
#include "FEmodel.h"

void FEElasticityMixedModel::GetSystem()
{
    int i,j,k,N=0,N1=0,N2=0;
    double temp;
    for(i=0;i<ms.Element_Num;i++)
    {
        j=dof1.Total_Num_PerElement[i]; k=dof2.Total_Num_PerElement[i];
        N+=j*j+2*j*k+3*k*k; N1+=j*j; N2+=k*k;
    }
    std::vector< Eigen::Triplet<double> > tripletList; //for constructing sparsmatrix
    tripletList.reserve(N);
    std::vector< Eigen::Triplet<double> > tripletList1; //for constructing sparsmatrix AL2
    tripletList1.reserve(N1);
    std::vector< Eigen::Triplet<double> > tripletList2; //for constructing sparsmatrix CH1
    tripletList2.reserve(N2);
    std::vector< Eigen::Triplet<double> > tripletListL2; //for constructing sparsmatrix CL1
    tripletList2.reserve(N2);
    double **LocAL2=new double*[polydim1];double **LocAtr=new double*[polydim1];double **LocAdiv=new double*[polydim1];
    double **LocB=new double*[polydim1];
    double **LocC1=new double*[polydim2];double **LocC2=new double*[polydim2];double **LocCL2=new double*[polydim2];
    double *LocF1=new double[polydim1];double *LocF2=new double[polydim2];
    for(j=0;j<polydim1;j++)
    {
        LocAL2[j]=new double[polydim1]; LocAtr[j]=new double[polydim1];   LocAdiv[j]=new double[polydim1];
        LocB[j]=new double[polydim2];
    }
    for(j=0;j<polydim2;j++)
    {
        LocC1[j]=new double[polydim2];LocC2[j]=new double[polydim2];LocCL2[j]=new double[polydim2];
    }
    for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
    {
        
        FE.GetMixedAL2(i, LocAL2); FE.GetMixedAL2trace(i, LocAtr); FE.GetMixedAdiv(i, LocAdiv);
        FE.GetMixedB(i, LocB);
        FE.GetMixedC1(i, LocC1);FE.GetMixedCH1(i, LocC2);FE.GetMixedCL2(i, LocCL2);
        FE.GetRHS1(i, pde.Source, LocF1);
        FE.GetRHS2(i, pde.Source, LocF2);
        for(j=0;j<polydim1;j++)//组装 A,B
        {
            RHS[ dof1.TotalD[i][j] ]+=gamma1*LocF1[j];//组装右端项F1
            for(k=0;k<polydim1;k++)
            {
                temp=1./(2*pde.mu)*(LocAL2[j][k]-pde.lambda/(2*pde.mu+2*pde.lambda)*LocAtr[j][k])+gamma1*LocAdiv[j][k];
                tripletList.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.TotalD[i][k],temp));
                tripletList1.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.TotalD[i][k],LocAL2[j][k]));
            }
            for(k=0;k<polydim2;k++)
            {
                tripletList.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.Dof_Num+dof2.TotalD[i][k],LocB[j][k]));
                tripletList.push_back(Eigen::Triplet<double>(dof1.Dof_Num+dof2.TotalD[i][k],dof1.TotalD[i][j],-LocB[j][k]));
            }
        }
        
        for(j=0;j<polydim2;j++)//组装 C1 F2
        {
            RHS[ dof1.Dof_Num+dof2.TotalD[i][j] ]+=LocF2[j];//组装右端项F1
            for(k=0;k<polydim2;k++)
            {
                tripletList.push_back(Eigen::Triplet<double>(dof1.Dof_Num+dof2.TotalD[i][j],dof1.Dof_Num+dof2.TotalD[i][k],gamma2*LocC1[j][k]));
                tripletList2.push_back(Eigen::Triplet<double>(dof2.TotalD[i][j],dof2.TotalD[i][k],LocC2[j][k]));
                tripletListL2.push_back(Eigen::Triplet<double>(dof2.TotalD[i][j],dof2.TotalD[i][k],LocCL2[j][k]));
            }
        }
    }
    for(i=0;i<ms.Edge_Num;i++)//组装刚度矩阵 C2
    {
        int E1ID=ms.EdgeElement[i][0],E2ID=ms.EdgeElement[i][1];
        if(E2ID==-1)
            continue;
        FE.GetMixedC2(i, LocC2);
        for(j=0;j<polydim2;j++)//组装 C2
            for(k=0;k<polydim2;k++)
            {
                tripletList.push_back(Eigen::Triplet<double>(dof1.Dof_Num+dof2.TotalD[E1ID][j],dof1.Dof_Num+dof2.TotalD[E2ID][k],gamma2*LocC2[j][k]));
                tripletList.push_back(Eigen::Triplet<double>(dof1.Dof_Num+dof2.TotalD[E2ID][k],dof1.Dof_Num+dof2.TotalD[E1ID][j],gamma2*LocC2[j][k]));
            }
    }
    StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
    AL2.setFromTriplets(tripletList1.begin(),tripletList1.end());
    CH1.setFromTriplets(tripletList2.begin(),tripletList2.end());
    CL2.setFromTriplets(tripletListL2.begin(),tripletListL2.end());
    
    for(j=0;j<polydim1;j++)
    {
        delete[] LocAL2[j]; delete[] LocAtr[j];   delete[] LocAdiv[j];  delete[] LocB[j];
    }
    for(j=0;j<polydim2;j++)
    {
       delete[] LocC1[j];delete[] LocC2[j];delete[] LocCL2[j];
    }
    delete[] LocAL2; delete[] LocAtr;  delete[] LocAdiv;  delete[] LocB;  delete[] LocC1;  delete[] LocC2; delete[] LocCL2;
}

void FEElasticityMixedModel::Solve()
{
    GetSystem();
    //deal with boundary condition for u
    int i,k,dofNum=dof1.Dof_Num+dof2.Dof_Num,BdofNum=ms.B_Edge_Num*dof2.Num_PerEdge,
    IdofNum=dofNum-BdofNum;;
    int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dofNum];
    double *BdofVal=new double[BdofNum];
    for(i=0;i<BdofNum;i++)
        BdofVal[i]=0;
    FE.GetBdof_BdofVal(pde.DirichletBoundary,Bdof,BdofVal);

    long N=StiffMatrix.nonZeros();
    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(N);
    
    VectorXd x0(dof1.Dof_Num+dof2.Dof_Num);
    x0.setZero();
    for(i=0;i<BdofNum;i++)
    {
        uh[Bdof[i]-dof1.Dof_Num]=BdofVal[i]; x0(Bdof[i])=BdofVal[i];
    }
    RHS-=StiffMatrix*x0;

    VectorXd b(IdofNum),x(IdofNum);
    int col=0,row;
    for(i=0;i<dofNum;i++) //get the relation between dof and inner dof
        if(i>=dof1.Dof_Num)
        {
            IdofID[i]=-1;
            if(dof2.DType[i-dof1.Dof_Num]==0)
            {
                b[col]=RHS[i];  Idof[col]=i; IdofID[i]=col;  col++;
            }
        }
        else
        {
            b[col]=RHS[i];  Idof[col]=i; IdofID[i]=col;  col++;
        }
    
    for(col=0;col<IdofNum;col++)
    {
        i=Idof[col];
        for(SparseMatrix<double>::InnerIterator it(StiffMatrix,i);it;++it)
        {
            k=(int) it.row();
            if(k>=dof1.Dof_Num)
            {
                if(dof2.DType[k-dof1.Dof_Num]==0)
                {
                    row=IdofID[k];
                    tripletList.push_back(Eigen::Triplet<double>(row,col,it.value()));
                }
            }
            else
            {
                row=IdofID[k];
                tripletList.push_back(Eigen::Triplet<double>(row,col,it.value()));
            }
        }
    }
    SparseMatrix<double> A(IdofNum,IdofNum);
    A.setFromTriplets(tripletList.begin(),tripletList.end());
    
    //direct method
    //Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    
    // iterative method
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    //Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    x=solver.solve(b);
    
    for(i=0;i<IdofNum;i++)
    {
        k=Idof[i];
        if(k<dof1.Dof_Num)
            sigmah[k]=x[i];
        else
            uh[k-dof1.Dof_Num]=x[i];
    }

    delete []Bdof; delete []BdofVal; delete []Idof; delete []IdofID;
}

double FEElasticityMixedModel::uH1Norm()
{
    double value=uI.transpose()*CH1*uI;
    return sqrt(value);
/*    double value=0;
    
   for(int ElemID=0;ElemID<ms.Element_Num;ElemID++)
    {
        int i,k;
    
    double uvalue[6],x[3],y[3];
    double *xpos=new double[4]; double *ypos=new double[4];//存储单元位置
    double xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<4;i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];   ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    for(i=0;i<4;i++)
    {
        x[0]=xE;x[1]=xpos[i];x[2]=xpos[(i+1)%4];
        y[0]=yE;y[1]=ypos[i];y[2]=ypos[(i+1)%4];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];     //积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            pde.Solution(jx,jy,uvalue);
            value+=(uvalue[1]*uvalue[1]+uvalue[2]*uvalue[2]+uvalue[4]*uvalue[4]+uvalue[5]*uvalue[5])*vol_K*jw;
        }
    }
    delete []xpos; delete []ypos;
    }
  return sqrt(value);
*/
}
double FEElasticityMixedModel::uhH1Norm()
{
    double value=uh.transpose()*CH1*uh;
    return sqrt(value);
}
double FEElasticityMixedModel::uhL2Norm()
{
    double value=uh.transpose()*CL2*uh;
    return sqrt(value);
}

double FEElasticityMixedModel::sigmaL2Norm()
{
    double value=sigmaI.transpose()*AL2*sigmah;
    return sqrt(value);
}

double FEElasticityMixedModel::sigmahL2Norm()
{
    double value=sigmah.transpose()*AL2*sigmah;
    return sqrt(value);
}

double FEElasticityMixedModel::uhH1Error()
{
    VectorXd X=uI-uh;
    double value=X.transpose()*CH1*X;
    return sqrt(value);
}

double FEElasticityMixedModel::uhL2Error()
{
    VectorXd X=uI-uh;
    double value=X.transpose()*CL2*X;
    return sqrt(value);
}

double FEElasticityMixedModel::uhMaxError()
{
    double value=0;
    int i;
    for(i=0;i<dof2.Dof_Num;i++)
        if(value<fabs(uI[i]-uh[i]))
            value=fabs(uI[i]-uh[i]);
    return value;
}

double FEElasticityMixedModel::sigmahL2Error()
{
    VectorXd X=sigmaI-sigmah;
    double value=X.transpose()*AL2*X;
    return sqrt(value);
}
double FEElasticityMixedModel::sigmahMaxError()
{
    double value=0;
    int i;
    for(i=0;i<dof1.Dof_Num;i++)
        if(value<fabs(sigmaI[i]-sigmah[i]))
            value=fabs(sigmaI[i]-sigmah[i]);
    return value;
}


