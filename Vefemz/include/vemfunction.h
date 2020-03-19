#ifndef vemfunction_h
#define vemfunction_h

#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::VectorXd;

#include "virtualelement.h"
#include "mathfunction.h"
#include "polynomialspace.h"

class VEMFunction
{
public:
	VEMFunction(VirtualElement &ve):SMS(ve.SMS),TQ(ve.TQ),ms(ve.ms),dof(ve.dof),VE(ve),pVector(ve.dof.Dof_Num)
	{
		size=dof.Dof_Num;
        pVector.setZero();
	}
    ~VEMFunction(){};
	void Interpolate(FunctionP u);
	void Interpolate(FunctionPt u,double t); //interpolation at t
	const double &operator[](int n) const;     //对[]进行重载，以返回第n个元素,不能改变数据
	double &operator[](int n);                 //对[]进行重载，以返回第n个元素，可以改变数据
	VEMFunction &operator=(const VEMFunction &pVer); //对=进行重载，以实现赋值，大小为零的向量可以被任意大小的向量赋值
    double GetL2ProjValue(double x,double y);//get L2 projection value at (x,y);

    double EnergyNorm(SparseMatrix<double> &A);
	double GetAverageValue();
	

    PolyMesh &ms;
    DegreeofFreedom &dof;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
	VirtualElement &VE;
    VectorXd pVector;
private:
    int size;

};




#endif
