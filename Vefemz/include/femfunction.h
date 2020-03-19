//
//  femfunction.h
//  Vefemz
//
//  Created by 张蓓 on 2019/11/10.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//

#ifndef femfunction_h
#define femfunction_h

#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::VectorXd;

#include "finiteelement.h"
#include "mathfunction.h"
#include "polynomialspace.h"

class FEMFunction
{
public:
    FEMFunction(FiniteElement &fe):SMS(fe.SMS),TQ(fe.TQ),ms(fe.ms),dof(fe.dof),FE(fe),pVector(fe.dof.Dof_Num)
    {
        size=dof.Dof_Num;
        pVector.setZero();
    }
    
    ~FEMFunction(){};
    void Interpolate(FunctionP u);
    void Interpolate(FunctionPt u,double t); //interpolation at t
    const double &operator[](int n) const;     //对[]进行重载，以返回第n个元素,不能改变数据
    double &operator[](int n);                 //对[]进行重载，以返回第n个元素，可以改变数据
    FEMFunction &operator=(const FEMFunction &pVer); //对=进行重载，以实现赋值，大小为零的向量可以被任意大小的向量赋值
    double GetValue(double x,double y);//get value at (x,y);
    
    double EnergyNorm(SparseMatrix<double> &A);
    double GetAverageValue();

    
    PolyMesh &ms;
    DegreeofFreedom &dof;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
    FiniteElement &FE;
    VectorXd pVector;
private:
    int size;
    
};
#endif /* femfunction_h */
