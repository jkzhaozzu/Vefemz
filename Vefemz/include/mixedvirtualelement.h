#ifndef mixedvirtualelement_h
#define mixedvirtualelement_h

#include "mesh.h"
#include "dof.h"
#include "quadrature.h"
#include "polynomialspace.h"
#include "mathfunction.h"
#include "virtualelement.h"

//mixed virtual element base class
class MixedVirtualElement
{
public:
    MixedVirtualElement(VirtualElement &ve1,VirtualElement &ve2):VE1(ve1),VE2(ve2) {};
    virtual void GetMixedB(int ElemID,double **B){};//-int_K pdiv v dx
    virtual void GetMixedA(int ElemID,double **A){};//int_Kgrad u:grad v dx
   
    VirtualElement &VE1;
    VirtualElement &VE2;
};

// H1-nonconforming p-order Stokes element with polynomial divergence
class MVEPkNC: public MixedVirtualElement
{
public:
    MVEPkNC();
    MVEPkNC(VirtualElement &ve1,VirtualElement &ve2): MixedVirtualElement(ve1,ve2) {;};
    
    ~MVEPkNC(){};
    virtual void GetMixedB(int ElemID,double **B);//-int_K pdiv v dx
    virtual void GetMixedA(int ElemID,double **A);//int_Kgrad u:grad v dx
    
};


// H1-nonconforming but H(div)-conforming p-order Stokes element with polynomial divergence
class MVEPkHdivNC: public MixedVirtualElement
{
public:
    MVEPkHdivNC();
    MVEPkHdivNC(VirtualElement &ve1,VirtualElement &ve2): MixedVirtualElement(ve1,ve2) {;};
    ~MVEPkHdivNC(){};
    void GetMixedB(int ElemID,double **B);//-int_K pdiv v dx
    void GetMixedA(int ElemID,double **A);//int_Kgrad u:grad v dx
    
};


#endif
