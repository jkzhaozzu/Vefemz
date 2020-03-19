#include <iostream>
using namespace std;
#include <iomanip>
#include "stdlib.h"
#include "stdio.h"
#include "time.h"

#include "mesh.h"
#include "dof.h"
#include "polynomialspace.h"
#include "quadrature.h"
#include "virtualelement.h"
#include "problemmodel.h"

void BiharmVEMorleyType()
{
	for(int meshID=1;meshID<6;meshID++)
	{
		
		clock_t start,end;
        double t1=0,t2=0;
		int k=2;//k>=2
		FILE *fp;
//		fp=GetMeshPointer(trimesh,meshID);
//		fp=GetMeshPointer(rectmesh,meshID);
		fp=GetMeshPointer(unipolymesh,meshID);
		PolyMesh ms(fp);		
		fclose(fp);	
		DegreeofFreedom dof(ms);
		ScaledMonomialSpace sms(k);
		TriangleQuadrature TQ(k+2);		
		VEMorleyType VE(k,sms,TQ,ms,dof);

/*		int ElemID=0;
		int i,j,m,NdofE=dof.Total_Num_PerElement[ElemID], pdim=(k+2)*(k+1)/2;
		double ** GH1=new double*[pdim];
		double ** BDH1=new double*[pdim];
		double ** BH1=new double*[pdim];
		double ** D=new double*[NdofE];
		double ** AH1=new double*[NdofE];
		for(i=0;i<pdim;i++)
		{
			GH1[i]=new double[pdim];
			BDH1[i]=new double[pdim];
			BH1[i]=new double[NdofE];
		}
		for(i=0;i<NdofE;i++)
		{
			D[i]=new double[pdim];
			AH1[i]=new double[NdofE];
		}
				
		VE.GetD(ElemID,D);VE.GetG_L2(ElemID,GH1);
		VE.GetB_L2(ElemID,BH1);VE.GetA_L2(ElemID,AH1);

		cout<<fixed << setprecision(7)<<endl;;
		for(i=0;i<NdofE;i++)
		{
			for(j=0;j<NdofE;j++)
			{
				cout<<fixed << setprecision(5)<<AH1[i][j]<<" ";
				BDH1[i][j]=0;
				for(m=0;m<NdofE;m++)
					BDH1[i][j]+=BH1[i][m]*D[m][j];
				cout<<"("<<i<<","<<j<<")"<<GH1[i][j]-BDH1[i][j]<<endl;

			}
			cout<<endl;
		}
*/
		cout<<fixed << setprecision(5);
		BiharmonicModel biharmonicequation(VE);
		
		start=clock();
		biharmonicequation.Solve();
		double uhEnergy=biharmonicequation.uh.EnergyNorm(biharmonicequation.StiffMatrix);
		double EnergyError=biharmonicequation.EnergyError();
		end=clock();

		t1=double(end-start)/CLOCKS_PER_SEC;
		cout<< "meshID= "<<meshID<<endl;
		cout<<"solving time: " << t1 <<" s"<< endl;
		cout<<uhEnergy<<"  "<<EnergyError<<"  "<<EnergyError/uhEnergy<<endl<<endl;

//		for(int i=0;i<dof.Dof_Num;i++)
//			cout<<"i="<<i<<" "<<biharmonicequation.uh[i]<<"  "<<biharmonicequation.RHS[i]<<endl;

	}
}


