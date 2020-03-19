#ifndef sparsesolve_h
#define sparsesolve_h

#include "stdlib.h"
#include "time.h"
#include <Eigen/Sparse>
#include "matrix.h"
#include "mathfunction.h"
#include "vemfunction.h"

void toEigenSparseSolve(Matrix &StiffMatrix,Vector &RHS,VEMFunction &uh);
void toEigenSparseSolve(Matrix &StiffMatrix,Vector &RHS,VEMFunction &uh,VEMFunction &ph);
void toEigenSparseSolveL02(Matrix &StiffMatrix,Vector &RHS,VEMFunction &uh,VEMFunction &ph);
#endif
