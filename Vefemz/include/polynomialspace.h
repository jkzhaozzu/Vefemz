//
//  polynomialspace.h
//  VEM2DforMac
//
//  Created by 张蓓 on 2019/1/3.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//

#ifndef polynomialspace_h
#define polynomialspace_h

//polynomial space base class
class PolynomialSpace
{
public:
	PolynomialSpace(int k):p(k){};
	~PolynomialSpace(){};
	
    
	virtual double GetValue(int i,double x,double y,double hE,double xE,double yE){return 0;};//get the i-th basis function value at (x,y)
	virtual double GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE){return 0;};//get i-th basis function (m,n)-derivative at (x,y)
	virtual void GetValue(int i,double x,double y,double hE,double xE,double yE,double *value){};
	virtual void GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE,double *value){};

    int p;    
};

class ScaledMonomialSpace: public PolynomialSpace
{
public:
	ScaledMonomialSpace(int k):PolynomialSpace(k)
	{
		int m=(p+1)*(p+2)/2;
		value=new double[m];
	};
    ~ScaledMonomialSpace()
	{
		delete [] value;
	};
    
    double GetValue(int i,double x,double y,double hE,double xE,double yE);//get the i-th basis function value at (x,y)
	double GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE);//get i-th basis function (m,n)-derivative at (x,y) 

private:
	double *value;
};

class ScaledMonomialSpaceV: public PolynomialSpace
{
public:
    ScaledMonomialSpaceV(int k):PolynomialSpace(k),SMS(k){};
			;
    ~ScaledMonomialSpaceV(){};
    
    void GetValue(int i,double x,double y,double hE,double xE,double yE,double *value);//get the i-th basis function value at (x,y)
	void GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE,double *value);//get i-th basis function (m,n)-derivative at (x,y)

private:
	ScaledMonomialSpace SMS;
    
};

// grad polynomial space
class GradPolynomialSpace: public PolynomialSpace
{
public:
    GradPolynomialSpace(int k):PolynomialSpace(k),SMS(k+1){};
			;
    ~GradPolynomialSpace(){};
    
    void GetValue(int i,double x,double y,double hE,double xE,double yE,double *value);//get the i-th basis function value at (x,y)
	void GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE,double *value);//get i-th basis function (m,n)-derivative at (x,y)

private:
	ScaledMonomialSpace SMS;
    
};

// curl polynomial space
class GradPolynomialOrthogonalSpace: public PolynomialSpace
{
public:
    GradPolynomialOrthogonalSpace(int k):PolynomialSpace(k),SMS(k-1){};
			;
    ~GradPolynomialOrthogonalSpace(){};
    
    void GetValue(int i,double x,double y,double hE,double xE,double yE,double *value);//get the i-th basis function value at (x,y)
	void GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE,double *value);//get i-th basis function (m,n)-derivative at (x,y)

private:
	ScaledMonomialSpace SMS;
    
};

#endif /* polynomialspace_h */
