#ifndef mathfuntions_h
#define mathfuntions_h

const double PI=4.0*atan(1.0);
typedef void (*FunctionP)(double x,double y,double *value);//function pointer type
typedef void (*FunctionPt)(double x,double y,double t,double *value);//function pointer type with variable t
//阶乘
int factorial(int m);

//Gauss消元法
void GaussSolve(int n,double **AA,double *bb,double *X);

int Sign(double x);

#endif