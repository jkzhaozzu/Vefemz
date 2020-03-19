
#include "stdio.h"
#include "stdlib.h"
#include "vemfunction.h"


void VEMFunction::Interpolate(FunctionP u)
{
	double *X=new double[dof.Dof_Num];
	for(int i=0;i<dof.Dof_Num;i++)
		X[i]=0;
	VE.GetDofVal(u,X);
	
	for(int i=0;i<dof.Dof_Num;i++)
		pVector[i]=X[i];
	delete [] X;
}

void VEMFunction::Interpolate(FunctionPt u,double t)
{
	double *X=new double[dof.Dof_Num];
	for(int i=0;i<dof.Dof_Num;i++)
		X[i]=0;
	VE.GetDofVal(u,t,X);
	
	for(int i=0;i<dof.Dof_Num;i++)
		pVector[i]=X[i];
	delete [] X;
}

double VEMFunction::GetAverageValue()
{
	int i,j;
	double value=0,measure=0;
	for(i=0;i<ms.Element_Num;i++)
	{
		double *X=new double[dof.Total_Num_PerElement[i]];
		for(j=0;j<dof.Total_Num_PerElement[i];j++)
			X[j]=pVector[dof.TotalD[i][j]];
		value+=VE.GetIntegralValue(i,X);
		measure+=ms.ElementMeasure[i];
		delete[] X;
	}
	value=value/measure;
	return value;
}

const double &VEMFunction::operator [](int n) const
{
    if(n>=0&&n<size)
    {
        return pVector[n];
    }
    else
    {
        printf("ERROR:向量坐标溢出\n");
        exit(1);
    }
}


double &VEMFunction::operator [](int n)
{
    if(n>=0&&n<size)
    {
        return pVector[n];
    }
    else
    {
        printf("ERROR:向量坐标溢出\n");
        exit(1);
    }
}

VEMFunction &VEMFunction::operator=(const VEMFunction &pVec)
{
    if(size==pVec.size)
    {
        pVector=pVec.pVector;
        return *this;
    }
    else
    {
        printf("ERROR:向量维数不同,且维数不同时左端不是零向量\n");
        exit(1);
    }
}

double VEMFunction::GetL2ProjValue(double x, double y)
{
    int ElemID=ms.FindElemID(x, y);
    int i,j,dofID;
    int p=VE.p,polydim=(p+2)*(p+1)/2,EdofNum=dof.Total_Num_PerElement[ElemID];
    double **PiL2=new double*[polydim];
    for(i=0;i<polydim;i++)
        PiL2[i]=new double[EdofNum];
    VE.GetPistar_H1(ElemID,PiL2);
    double value=0,hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<EdofNum;i++)
    {
        dofID=dof.TotalD[ElemID][i];
        for(j=0;j<polydim;j++)
            value+=pVector[dofID]*PiL2[j][i]*SMS.GetValue(j, x, y, hE, xE, yE);
    }
    for(i=0;i<polydim;i++)
        delete []PiL2[i];
    delete []PiL2;
    
    return value;
}

double VEMFunction::EnergyNorm(Eigen::SparseMatrix<double> &A)
{
    double value=pVector.transpose()*A*pVector;
    return sqrt(value);
}


