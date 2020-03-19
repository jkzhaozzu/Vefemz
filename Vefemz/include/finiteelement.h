//
//  finiteelement.h
//  Vefemz
//
//  Created by 张蓓 on 2019/11/10.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//

#ifndef finiteelement_h
#define finiteelement_h

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "mesh.h"
#include "dof.h"
#include "quadrature.h"
#include "polynomialspace.h"
#include "mathfunction.h"

//finite element base class
class FiniteElement
{
public:
    FiniteElement(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):
    p(k),SMS(ps),TQ(tq),ms(mesh),dof(Dof) { };
    
    virtual void GetD(int ElemID,double **D){};
    virtual void GetBase(int ElemID,double **BF){};
    virtual void GetG_L2(int ElemID,double **G){};
    virtual void GetA_L2(int ElemID,double **A){};
    virtual void GetG_H1(int ElemID,double **G){};
    virtual void GetA_H1(int ElemID,double **A){};
    virtual void GetG_H2(int ElemID,double **G){};
    virtual void GetA_H2(int ElemID,double **A){};
    virtual void GetB_div(int ElemID,double **B){};
    virtual void GetG_div(int ElemID,double **G){};
    virtual void GetPistar_div(int ElemID,double **Pistar){};
    virtual void GetA_div(int ElemID,double **A){};
    virtual double GetIntegralValue(int ElemID,double *X){return 0;};
    virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval){};// get boundary dof ID and corresponding dof value
    virtual void GetDofVal(FunctionP u,double *uI){};
    virtual void GetRHS(int ElemID,FunctionP Source,double *LocF){};
    virtual void GetBdof_BdofVal(FunctionPt BFunc,double t,int *Bdof,double *Bdofval){};// get boundary dof ID and corresponding dof value at t
    virtual void GetDofVal(FunctionPt u,double t,double *uI){};
    virtual void GetRHS(int ElemID,FunctionPt Source,double t,double *LocF){};
    
    virtual void GetdxB_L2(int ElemID,double **B){};
    virtual void GetdyB_L2(int ElemID,double **B){};
    virtual void GetdxPistar_L2(int ElemID,double **Pistar){};
    virtual void GetdyPistar_L2(int ElemID,double **Pistar){};
    
    int p; //the degree of element
    PolyMesh &ms;
    DegreeofFreedom &dof;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
};

//C0 conforming p-order element
class FEPkC0: public FiniteElement
{
public:
    FEPkC0();
    FEPkC0(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):FiniteElement(k,ps,tq,mesh,Dof)
    {
        polydim=(p+2)*(p+1)/2;
        dof.Num_PerNode=1; dof.Num_PerEdge=p-1; dof.Num_PerElement=(p-1)*(p-2)/2;
        dof.GetDofInfo();  dof.MallocDofMemory();  dof.Construct();
    };
    ~FEPkC0(){};
    
    
    virtual void GetD(int ElemID,double **D);
    virtual void GetBase(int ElemID,double **BF);
    virtual void GetG_L2(int ElemID,double **G);
    virtual void GetA_L2(int ElemID,double **A);
    virtual void GetG_H1(int ElemID,double **G);
    virtual void GetA_H1(int ElemID,double **A);
    virtual void GetRHS(int ElemID,FunctionP Source,double *LocF);
    virtual void GetDofVal(FunctionP u,double *uI);
    virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
/*    virtual void GetRHS(int ElemID,FunctionPt Source,double t,double *LocF);
    virtual void GetDofVal(FunctionPt u,double t,double *uI);
    virtual void GetBdof_BdofVal(FunctionPt BFunc,double t,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value at t
*/
private:
    void GetRHSL2B(int ElemID,FunctionP Source,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID
//    void GetRHSL2B(int ElemID,FunctionPt Source,double t,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID at t
    int polydim;//the dim of p-order polynomial space
 
};

#endif /* finiteelement_h */

