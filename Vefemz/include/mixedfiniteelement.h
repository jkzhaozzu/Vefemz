//
//  mixedfiniteelement.h
//  Vefemz
//
//  Created by 张蓓 on 2020/2/17.
//  Copyright © 2020年 Jikun Zhao. All rights reserved.
//

#ifndef mixedfiniteelement_h
#define mixedfiniteelement_h

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "mesh.h"
#include "dof.h"
#include "quadrature.h"
#include "polynomialspace.h"
#include "mathfunction.h"

//stablized mixed FE for linear elasticity on rectangular meshes

class MFEP1NCRS
{
public:
    MFEP1NCRS();
    MFEP1NCRS(TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof1,DegreeofFreedom & Dof2):
    TQ(tq),ms(mesh),dof1(Dof1),dof2(Dof2)
    {
        polydim1=5;polydim2=8;
        dof1.Num_PerNode=0; dof1.Num_PerEdge=0; dof1.Num_PerElement=5;
        dof1.GetDofInfo();  dof1.MallocDofMemory();  dof1.Construct();
        dof2.Num_PerNode=0; dof2.Num_PerEdge=2; dof2.Num_PerElement=0;
        dof2.GetDofInfo();  dof2.MallocDofMemory();  dof2.Construct();
    };
    ~MFEP1NCRS(){};
    
    void GetMixedAL2(int ElemID,double **A);//int_K sigma:tau dx
    void GetMixedAL2trace(int ElemID,double **Atr);//int_K tr(sigma):tr(tau) dx
    void GetMixedAdiv(int ElemID,double **A);//h_K^2*int_Kdiv(sigma)div(tau) dx
    void GetMixedB(int ElemID,double **B);//-int_K tau:epsilon(v) dx
    void GetMixedC1(int ElemID,double **C);//h_E^-1\int_E uv ds
    void GetMixedC2(int EdgeID,double **C);//-h_E^-1\int_E u1v2 ds
   // void GetMixedC3(int EdgeID,double **C);//-h_E^-1\int_E u2v1 ds
    void GetMixedCH1(int EdgeID,double **C);//(grad u,grad v)_K
    void GetMixedCL2(int EdgeID,double **C);//(u,v)_K
    void GetRHS1(int ElemID,FunctionP Source,double *LocF);// h_K^2(f,div tau)_K
    void GetRHS2(int ElemID,FunctionP Source,double *LocF);// (f,vh)_K
    
    void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
    void GetDof1Val(FunctionP sigma,VectorXd &sigmaI); //get dof value of sigma
    void GetDof2Val(FunctionP u,VectorXd &uI); //get dof value of u
    
    PolyMesh &ms;
    DegreeofFreedom &dof1,&dof2;
    TriangleQuadrature &TQ;
    
private:
    double polydim1,polydim2;
    void GetBasisFuntion1(int ElemID,int i,double x,double y,double *value);//get basis function from first space
    void GetBasisFuntion2(int ElemID,int i,double x,double y,double *value);//get basis function from second space
    void CheckBasisFunc(int ElemID,double **C);//for u
    
};


#endif /* mixedfiniteelement_h */
