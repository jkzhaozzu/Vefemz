#include "mixedvirtualelement.h"

void MVEPkNC::GetMixedB(int ElemID,double **B)
{
	int i,j,k,EdofN1=VE1.dof.Total_Num_PerElement[ElemID],EdofN2=VE2.dof.Total_Num_PerElement[ElemID];
	double **Pistar=new double*[EdofN2];
	double **BT=new double*[EdofN2];
	for(i=0;i<EdofN2;i++)
	{
		BT[i]=new double[EdofN1];Pistar[i]=new double[EdofN2];
	}
	VE1.GetB_div(ElemID,BT); VE2.GetPistar_L2(ElemID,Pistar);
	for(i=0;i<EdofN1;i++)
		for(j=0;j<EdofN2;j++)
		{
			B[i][j]=0;
			for(k=0;k<EdofN2;k++)	
				B[i][j]-=Pistar[k][j]*BT[k][i];
		}

	for(i=0;i<EdofN2;i++)
	{
		delete []BT[i]; delete []Pistar[i];
	}
	delete []BT; delete []Pistar;
}

void MVEPkNC::GetMixedA(int ElemID,double **A)
{
	VE1.GetA_H1(ElemID,A);
}