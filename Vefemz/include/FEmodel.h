//
//  FEmodel.h
//  Vefemz
//
//  Created by 张蓓 on 2019/11/10.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//

#ifndef FEmodel_h
#define FEmodel_h

#include "finiteelement.h"
#include "mixedfiniteelement.h"
#include "femfunction.h"
#include "problemmodel.h"


// Poisson solver
class FEPoissonModel
{
public:
    FEPoissonModel(FiniteElement &fe):SMS(fe.SMS),TQ(fe.TQ),ms(fe.ms),dof(fe.dof),FE(fe),uh(fe),uI(fe),StiffMatrix(dof.Dof_Num,dof.Dof_Num),RHS(dof.Dof_Num)
    {
        polydim=(FE.p+2)*(FE.p+1)/2; RHS.setZero();
        uI.Interpolate(pde.Solution);
    }
    ~FEPoissonModel(){};
    void GetLocRHS(int ElemID,double *LocF);//get local RHS on element ElemID
    void GetSystem();// get StiffMatrix and RHS
    void Solve();//deal with BD and get uh

    double EnergyError();
    double MaxError();

    int polydim;
    PolyMesh &ms;
    DegreeofFreedom &dof;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
    FiniteElement &FE;
    FEMFunction uh;//numerical solution
    FEMFunction uI;//interpolation of solution
    PoissonProblem pde;
    SparseMatrix<double> StiffMatrix;
    VectorXd RHS;
    
};


// Elasticity Mixed FEM solver with stabilization
class FEElasticityMixedModel
{
public:
    FEElasticityMixedModel(MFEP1NCRS &fe):
    TQ(fe.TQ),ms(fe.ms),dof1(fe.dof1),dof2(fe.dof2),FE(fe),sigmah(dof1.Dof_Num),sigmaI(dof1.Dof_Num),uh(dof2.Dof_Num),uI(dof2.Dof_Num),
    StiffMatrix(dof1.Dof_Num+dof2.Dof_Num,dof1.Dof_Num+dof2.Dof_Num),AL2(dof1.Dof_Num,dof1.Dof_Num),
    CH1(dof2.Dof_Num,dof2.Dof_Num),CL2(dof2.Dof_Num,dof2.Dof_Num),
    RHS(dof1.Dof_Num+dof2.Dof_Num)
    {
        polydim1=5; polydim2=8; gamma1=0.05;gamma2=1;RHS.setZero();
        sigmah.setZero();sigmaI.setZero();uh.setZero();uI.setZero();
        FE.GetDof1Val(pde.Stress, sigmaI); FE.GetDof2Val(pde.Solution, uI);
    }
    ~FEElasticityMixedModel(){};
    void GetSystem();// get StiffMatrix and RHS
    void Solve();//deal with BD and get uh
    double uH1Norm();
    double uhH1Norm();
    double uhL2Norm();
    double sigmaL2Norm();
    double sigmahL2Norm();
    
    double uhH1Error();
    double uhL2Error();
    double uhMaxError();
    double sigmahL2Error();
    double sigmahMaxError();
    
    int polydim1,polydim2;
    double gamma1,gamma2;
    PolyMesh &ms;
    DegreeofFreedom &dof1,&dof2;
    TriangleQuadrature &TQ;
    MFEP1NCRS &FE;
    VectorXd sigmah,sigmaI,uh,uI;//uh:numerical solution, uI:interpolation of solution
    ElasticityProblem pde;
    SparseMatrix<double> StiffMatrix,AL2,CH1,CL2;
    VectorXd RHS;
};


#endif /* FEmodel_h */
