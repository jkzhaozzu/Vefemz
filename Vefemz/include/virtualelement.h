//
//  virtualelement.h
//  VEM2DforMac
//
//  Created by 张蓓 on 2019/1/7.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//

#ifndef virtualelement_h
#define virtualelement_h

#include "mesh.h"
#include "dof.h"
#include "quadrature.h"
#include "polynomialspace.h"
#include "mathfunction.h"


//virtual element base class
class VirtualElement
{
public:
	VirtualElement(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):
    p(k),SMS(ps),TQ(tq),ms(mesh),dof(Dof) {	};

	virtual void GetB_H1(int ElemID,double **B){};
	virtual void GetD(int ElemID,double **D){};
	virtual void GetG_H1(int ElemID,double **G){};
	virtual void GetPistar_H1(int ElemID,double **Pistar){};
	virtual void GetPi_H1(int ElemID,double **Pi){};
	virtual void GetA_H1(int ElemID,double **A){};
	virtual void GetB_H2(int ElemID,double **B){};
	virtual void GetG_H2(int ElemID,double **G){};
	virtual void GetPistar_H2(int ElemID,double **Pistar){};
	virtual void GetPi_H2(int ElemID,double **Pi){};
	virtual void GetA_H2(int ElemID,double **A){};
	virtual void GetB_L2(int ElemID,double **B){};
	virtual void GetG_L2(int ElemID,double **G){};
	virtual void GetPistar_L2(int ElemID,double **Pistar){};
	virtual void GetPi_L2(int ElemID,double **Pi){};	
	virtual void GetA_L2(int ElemID,double **A){};
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
class VEPkC0: public VirtualElement
{
public:
    VEPkC0();
    VEPkC0(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):VirtualElement(k,ps,tq,mesh,Dof)
    {
		polydim=(p+2)*(p+1)/2;
		dof.Num_PerNode=1; dof.Num_PerEdge=p-1; dof.Num_PerElement=p*(p-1)/2;
        dof.GetDofInfo();  dof.MallocDofMemory();  dof.Construct();	
    };
    ~VEPkC0(){};
    
    virtual void GetB_H1(int ElemID,double **B);
    virtual void GetD(int ElemID,double **D);
    virtual void GetG_H1(int ElemID,double **G);
	virtual void GetPistar_H1(int ElemID,double **Pistar);
	virtual void GetPi_H1(int ElemID,double **Pi);
	virtual void GetB_L2(int ElemID,double **B);
	virtual void GetG_L2(int ElemID,double **G);
	virtual void GetPistar_L2(int ElemID,double **Pistar);
	virtual void GetPi_L2(int ElemID,double **Pi);
	virtual void GetA_H1(int ElemID,double **A);
	virtual void GetA_L2(int ElemID,double **A);
	virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
	virtual void GetDofVal(FunctionP u,double *uI);
	virtual void GetRHS(int ElemID,FunctionP Source,double *LocF);
	virtual void GetBdof_BdofVal(FunctionPt BFunc,double t,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value at t
	virtual void GetDofVal(FunctionPt u,double t,double *uI);
	virtual void GetRHS(int ElemID,FunctionPt Source,double t,double *LocF);

private:
	void GetRHSL2B(int ElemID,FunctionP Source,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID
	void GetRHSL2B(int ElemID,FunctionPt Source,double t,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID at t
	int polydim;//the dim of p-order polynomial space
	
  
};


//H1-nonconforming p-order element 
class VEPkNC: public VirtualElement
{
public:
    VEPkNC();
    VEPkNC(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof): VirtualElement(k,ps,tq,mesh,Dof)
    {
		polydim=(p+2)*(p+1)/2;
		dof.Num_PerNode=0; dof.Num_PerEdge=p;  dof.Num_PerElement=p*(p-1)/2;
        dof.GetDofInfo(); dof.MallocDofMemory(); dof.Construct();
    };
    ~VEPkNC(){};
    
    virtual void GetB_H1(int ElemID,double **B);
    virtual void GetD(int ElemID,double **D);
    virtual void GetG_H1(int ElemID,double **G);
	virtual void GetPistar_H1(int ElemID,double **Pistar);
	virtual void GetPi_H1(int ElemID,double **Pi);
	virtual void GetB_L2(int ElemID,double **B);
	virtual void GetG_L2(int ElemID,double **G);
	virtual void GetPistar_L2(int ElemID,double **Pistar);
	virtual void GetPi_L2(int ElemID,double **Pi);
	virtual void GetA_H1(int ElemID,double **A);
	virtual void GetA_L2(int ElemID,double **A);
	virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
	virtual void GetDofVal(FunctionP u,double *uI);
	virtual void GetRHS(int ElemID,FunctionP Source,double *LocF);
	virtual void GetBdof_BdofVal(FunctionPt BFunc,double t,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value at t
	virtual void GetDofVal(FunctionPt u,double t,double *uI);
	virtual void GetRHS(int ElemID,FunctionPt Source,double t,double *LocF);

private:
	void GetRHSL2B(int ElemID,FunctionP Source,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID
	void GetRHSL2B(int ElemID,FunctionPt Source,double t,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID at t
	int polydim;//the dim of p-order polynomial space
    
};

//H1-nonconforming p-order element with node dof
class VEPkNCNodeType: public VirtualElement
{
public:
    VEPkNCNodeType();
    VEPkNCNodeType(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof): VirtualElement(k,ps,tq,mesh,Dof)
    {
		polydim=(p+2)*(p+1)/2;
		dof.Num_PerNode=0; dof.Num_PerEdge=p;  dof.Num_PerElement=p*(p-1)/2;
        dof.GetDofInfo();  dof.MallocDofMemory(); dof.Construct();
    };
    ~VEPkNCNodeType(){};
    
    virtual void GetB_H1(int ElemID,double **B);
    virtual void GetD(int ElemID,double **D);
    virtual void GetG_H1(int ElemID,double **G);
	virtual void GetPistar_H1(int ElemID,double **Pistar);
	virtual void GetPi_H1(int ElemID,double **Pi);
	virtual void GetB_L2(int ElemID,double **B);
	virtual void GetG_L2(int ElemID,double **G);
	virtual void GetPistar_L2(int ElemID,double **Pistar);
	virtual void GetPi_L2(int ElemID,double **Pi);
	virtual void GetA_H1(int ElemID,double **A);
	virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
	virtual void GetDofVal(FunctionP u,double *uI);

private:
	int polydim;//the dim of p-order polynomial space
    
};

// H1-nonconforming p-order vector-valued element with polynomial divergence
class VEPkNCV:public VirtualElement
{
public:
    VEPkNCV();
	VEPkNCV(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):VirtualElement(k,ps,tq,mesh,Dof)
    {
		polydim=(p+2)*(p+1);
        dof.Num_PerNode=0; dof.Num_PerEdge=2*p; dof.Num_PerElement=p*(p-1);         
        dof.GetDofInfo(); dof.MallocDofMemory(); dof.Construct();	
    };
    ~VEPkNCV(){};
    
	virtual void GetB_H1(int ElemID,double **B);
    virtual void GetD(int ElemID,double **D);
    virtual void GetG_H1(int ElemID,double **G);
	virtual void GetPistar_H1(int ElemID,double **Pistar);
	virtual void GetPi_H1(int ElemID,double **Pi);
	virtual void GetB_L2(int ElemID,double **B);
	virtual void GetG_L2(int ElemID,double **G);
	virtual void GetPistar_L2(int ElemID,double **Pistar);
	virtual void GetPi_L2(int ElemID,double **Pi);
	virtual void GetA_H1(int ElemID,double **A);
	virtual void GetA_L2(int ElemID,double **A);
	virtual void GetB_div(int ElemID,double **B);
	virtual void GetG_div(int ElemID,double **G);
	virtual void GetPistar_div(int ElemID,double **Pistar);
	virtual void GetA_div(int ElemID,double **A);
	
	virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
	virtual void GetDofVal(FunctionP u,double *uI);
	virtual void GetRHS(int ElemID,FunctionP Source,double *LocF);

	virtual void GetdxB_L2(int ElemID,double **B);
	virtual void GetdyB_L2(int ElemID,double **B);
	virtual void GetdxPistar_L2(int ElemID,double **Pistar);
	virtual void GetdyPistar_L2(int ElemID,double **Pistar);
private:
	void GetRHSL2B(int ElemID,FunctionP Source,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID
  
	void L2ProjecttoEdgePoly(int ElemID,int EdgeID,double **Pistar);//projection of v\cdot n on edge EdgeID of Element ElemID
	void PolynomialDecomposition(int ElemID,double **PistarGrad,double **PistarGradOrth);

	int polydim;//the dim of p-order polynomial vector-valued space
};

// H1-nonconforming but H(div)-conforming p-order vector-valued element with polynomial divergence
class VEPkHdivNCV:public VirtualElement
{
public:
    VEPkHdivNCV();
    VEPkHdivNCV(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):VirtualElement(k,ps,tq,mesh,Dof)
    {
        polydim=(p+2)*(p+1);
        dof.Num_PerNode=0; dof.Num_PerEdge=2*p+1; dof.Num_PerElement=p*(p-1);
        dof.GetDofInfo(); dof.MallocDofMemory(); dof.Construct();
    };
    ~VEPkHdivNCV(){};
    
    virtual void GetB_H1(int ElemID,double **B);
    virtual void GetD(int ElemID,double **D);
    virtual void GetG_H1(int ElemID,double **G);
    virtual void GetPistar_H1(int ElemID,double **Pistar);
    virtual void GetPi_H1(int ElemID,double **Pi);
    virtual void GetB_L2(int ElemID,double **B);
    virtual void GetG_L2(int ElemID,double **G);
    
    virtual void GetPistar_L2(int ElemID,double **Pistar);
    virtual void GetPi_L2(int ElemID,double **Pi);
    virtual void GetA_H1(int ElemID,double **A);
    virtual void GetA_L2(int ElemID,double **A);
    
    virtual void GetB_div(int ElemID,double **B);
    virtual void GetG_div(int ElemID,double **G);
    virtual void GetPistar_div(int ElemID,double **Pistar);
    virtual void GetA_div(int ElemID,double **A);
    

    virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
    virtual void GetDofVal(FunctionP u,double *uI);
    virtual void GetRHS(int ElemID,FunctionP Source,double *LocF);
 /*
    virtual void GetdxB_L2(int ElemID,double **B);
    virtual void GetdyB_L2(int ElemID,double **B);
    virtual void GetdxPistar_L2(int ElemID,double **Pistar);
    virtual void GetdyPistar_L2(int ElemID,double **Pistar);
*/
private:
    void GetRHSL2B(int ElemID,FunctionP Source,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID
    void L2ProjecttoEdgePoly(int ElemID,int EdgeID,double **Pistar);//projection of v\cdot n on edge EdgeID of Element ElemID
    void PolynomialDecomposition(int ElemID,double **PistarGrad,double **PistarGradOrth);

    int polydim;//the dim of p-order polynomial vector-valued space
};

//DG p-order element 
class VEPkDG: public VirtualElement
{
public:
    VEPkDG();
    VEPkDG(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):VirtualElement(k,ps,tq,mesh,Dof)
    {
		polydim=(p+2)*(p+1)/2;
		dof.Num_PerNode=0; dof.Num_PerEdge=0; dof.Num_PerElement=(p+2)*(p+1)/2;
        dof.GetDofInfo();  dof.MallocDofMemory();  dof.Construct();	
    };
    ~VEPkDG(){};
    virtual void GetG_L2(int ElemID,double **G);
	virtual void GetB_L2(int ElemID,double **G);
	virtual void GetPistar_L2(int ElemID,double **Pistar);
	virtual void GetA_L2(int ElemID,double **A);
	virtual double GetIntegralValue(int ElemID,double *X);
//	virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
	virtual void GetDofVal(FunctionP u,double *uI);
//	virtual void GetRHS(int ElemID,FunctionP Source,double *LocF);

private:
	void GetRHSL2B(int ElemID,FunctionP Source,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID
	int polydim;//the dim of p-order polynomial space
	
  
};

//H2 conforming p-order element
class VEPkH2: public VirtualElement
{
public:
    VEPkH2();
    VEPkH2(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):VirtualElement(k,ps,tq,mesh,Dof)
    {
        polydim=(p+2)*(p+1)/2;
        int r=p;
        if(r<3) r=3;
        dof.Num_PerNode=3; dof.Num_PerEdge=r-3+p-2; dof.Num_PerElement=(p-2)*(p-3)/2;
        dof.GetDofInfo();  dof.MallocDofMemory();  dof.Construct();
        Getvcharlength();
    };
    ~VEPkH2(){ delete []vcharlength; };
    
    virtual void GetB_H2(int ElemID,double **B);
    virtual void GetD(int ElemID,double **D);
    virtual void GetG_H2(int ElemID,double **G);
    virtual void GetPistar_H2(int ElemID,double **Pistar);
    virtual void GetPi_H2(int ElemID,double **Pi);
    virtual void GetA_H2(int ElemID,double **A);
/*    virtual void GetB_H1(int ElemID,double **B);
    virtual void GetG_H1(int ElemID,double **G);
    virtual void GetPistar_H1(int ElemID,double **Pistar);
    virtual void GetPi_H1(int ElemID,double **Pi);
    virtual void GetA_H1(int ElemID,double **A);
*/
    virtual void GetG_L2(int ElemID,double **G);

    virtual void GetRHS(int ElemID,FunctionP Source,double *LocF);//(P_k-2f,v)
    virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
    virtual void GetDofVal(FunctionP u,double *uI);


private:
//    void L2ProjecttoEdgePoly(int ElemID,int EdgeID,double **Pistar);//projection of v on edge EdgeID of Element ElemID
    void GetRHSL2B(int ElemID,FunctionP Source,double *LocFB);//get the quadrature of RHS with space P_k-2 on element ElemID
    void Getvcharlength();
    int polydim;//the dim of p-order polynomial space
    double *vcharlength;//the element diameter average w.r.t. each vertex
    
    
};

//C0-continous H2 nonconforming p-order element 
class VEC0PkH2NC: public VirtualElement
{
public:
    VEC0PkH2NC();
    VEC0PkH2NC(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):VirtualElement(k,ps,tq,mesh,Dof)
    {
		polydim=(p+2)*(p+1)/2;
		dof.Num_PerNode=1; dof.Num_PerEdge=2*(p-1); dof.Num_PerElement=(p-2)*(p-3)/2;
        dof.GetDofInfo();  dof.MallocDofMemory();  dof.Construct();	
    };
    ~VEC0PkH2NC(){};
    
	virtual void GetB_H2(int ElemID,double **B);
	virtual void GetD(int ElemID,double **D);
    virtual void GetG_H2(int ElemID,double **G);
	virtual void GetPistar_H2(int ElemID,double **Pistar);
	virtual void GetPi_H2(int ElemID,double **Pi);
	virtual void GetA_H2(int ElemID,double **A);
	virtual void GetB_H1(int ElemID,double **B);
	virtual void GetG_H1(int ElemID,double **G);
	virtual void GetPistar_H1(int ElemID,double **Pistar);
	virtual void GetPi_H1(int ElemID,double **Pi);
	virtual void GetA_H1(int ElemID,double **A);
    
    virtual void GetB_L2(int ElemID,double **B);
    virtual void GetG_L2(int ElemID,double **G);
    virtual void GetPistar_L2(int ElemID,double **Pistar);
    virtual void GetPi_L2(int ElemID,double **Pi);
    virtual void GetA_L2(int ElemID,double **A);
	
	virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
	virtual void GetDofVal(FunctionP u,double *uI);
	virtual void GetRHS(int ElemID,FunctionP Source,double *LocF);//(P_k-2f,v)
    
    virtual void GetBdof_BdofVal(FunctionPt BFunc,double t,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value at t
    virtual void GetDofVal(FunctionPt u,double t,double *uI);
    virtual void GetRHS(int ElemID,FunctionPt Source,double t,double *LocF);

private:
	void L2ProjecttoEdgePoly(int ElemID,int EdgeID,double **Pistar);//projection of v on edge EdgeID of Element ElemID
	void GetRHSL2B(int ElemID,FunctionP Source,double *LocFB);//get the quadrature of RHS with space P_k-2 on element ElemID
    void GetRHSL2B(int ElemID,FunctionPt Source,double t,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID at t
	int polydim;//the dim of p-order polynomial space
	
  
};

//Morley type element 
class VEMorleyType: public VirtualElement
{
public:
    VEMorleyType();
    VEMorleyType(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):VirtualElement(k,ps,tq,mesh,Dof)
    {
		polydim=(p+2)*(p+1)/2;
		dof.Num_PerNode=1; dof.Num_PerEdge=2*p-3; dof.Num_PerElement=(p-2)*(p-3)/2;
        dof.GetDofInfo();  dof.MallocDofMemory();  dof.Construct();	
    };
    ~VEMorleyType(){};
    
	virtual void GetB_H2(int ElemID,double **B);
	virtual void GetD(int ElemID,double **D);
	virtual void GetG_H2(int ElemID,double **G);
	virtual void GetPistar_H2(int ElemID,double **Pistar);
	virtual void GetPi_H2(int ElemID,double **Pi);
	virtual void GetA_H2(int ElemID,double **A);
	
    virtual void GetB_L2(int ElemID,double **B);
    virtual void GetG_L2(int ElemID,double **G);
    virtual void GetPistar_L2(int ElemID,double **Pistar);
    virtual void GetPi_L2(int ElemID,double **Pi);
    virtual void GetA_L2(int ElemID,double **A);
/*	virtual void GetB_H1(int ElemID,double **B);
	virtual void GetG_H1(int ElemID,double **G);
	virtual void GetPistar_H1(int ElemID,double **Pistar);
	virtual void GetPi_H1(int ElemID,double **Pi);
	virtual void GetA_H1(int ElemID,double **A);
*/	
	virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
	virtual void GetDofVal(FunctionP u,double *uI);
	virtual void GetRHS(int ElemID,FunctionP Source,double *LocF);//(P_k-2f,v)
private:
//	void L2ProjecttoEdgePoly(int ElemID,int EdgeID,double **Pistar);//projection of v on edge EdgeID of Element ElemID
	void GetRHSL2B(int ElemID,FunctionP Source,double *LocFB);//get the quadrature of RHS with space P_k-2 on element ElemID
	int polydim;//the dim of p-order polynomial space
	
  
};

// H1-nonconforming p-order vector-valued element with non-polynomial divergence
class VEPkNCVNonPolyDiv:public VirtualElement
{
public:
    VEPkNCVNonPolyDiv();
	VEPkNCVNonPolyDiv(int k,PolynomialSpace & ps, TriangleQuadrature &tq, PolyMesh & mesh,DegreeofFreedom & Dof):
		VirtualElement(k,ps,tq,mesh,Dof)
    {
		polydim=(p+2)*(p+1);
		dof.Num_PerNode=0; dof.Num_PerEdge=2*p; dof.Num_PerElement=p*(p-1);
        dof.GetDofInfo();  dof.MallocDofMemory(); dof.Construct();	
    };
    ~VEPkNCVNonPolyDiv(){};
    
	virtual void GetB_H1(int ElemID,double **B);
    virtual void GetD(int ElemID,double **D);
    virtual void GetG_H1(int ElemID,double **G);
	virtual void GetPistar_H1(int ElemID,double **Pistar);
	virtual void GetPi_H1(int ElemID,double **Pi);
	virtual void GetB_L2(int ElemID,double **B);
	virtual void GetG_L2(int ElemID,double **G);
	virtual void GetPistar_L2(int ElemID,double **Pistar);
	virtual void GetPi_L2(int ElemID,double **Pi);
	virtual void GetA_H1(int ElemID,double **A);
	virtual void GetB_div(int ElemID,double **B);
	virtual void GetG_div(int ElemID,double **G);
	virtual void GetPistar_div(int ElemID,double **Pistar);
	virtual void GetA_div(int ElemID,double **A);
/*	virtual void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value
	virtual void GetDofVal(FunctionP u,double *uI);
 */ 
private:
	int polydim;//the dim of p-order polynomial vector-valued space
};


#endif /* virtualelement_h */
