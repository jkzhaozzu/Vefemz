#include <iostream>
using namespace std;
#include "stdlib.h"
#include "stdio.h"
#include "time.h"

#include "mesh.h"
#include "dof.h"
#include "polynomialspace.h"
#include "quadrature.h"
#include "virtualelement.h"
#include "problemmodel.h"



void ParabolicVEPkC0()
{

	for(int meshID=1;meshID<8;meshID++)
	{
		
		clock_t start,end;
        double t1=0;
		int k=1;
		FILE *fp;
//		fp=GetMeshPointer(trimesh,meshID);
//		fp=GetMeshPointer(rectmesh,meshID);
		fp=GetMeshPointer(unipolymesh,meshID);
		PolyMesh ms(fp);		
		fclose(fp);	
		DegreeofFreedom dof(ms);
		ScaledMonomialSpace sms(k);
		TriangleQuadrature TQ(k+2);		
		VEPkC0 VE(k,sms,TQ,ms,dof);
		int dtn=100;//iteration number
		ParabolicModel parabolic(VE,dtn);
		
		start=clock();
		parabolic.Solve();
		double uhH1=parabolic.uh.EnergyNorm(parabolic.MatrixH1),uhL2=parabolic.uh.EnergyNorm(parabolic.MatrixL2);
		double ErrorH1=parabolic.ErrorH1(),ErrorL2=parabolic.ErrorL2();
		end=clock();

		t1=double(end-start)/CLOCKS_PER_SEC;
		cout<< "meshID= "<<meshID<<endl;
		cout<<"solving time: " << t1 <<" s"<< endl;
		cout<<"H1error:"<<uhH1<<"  "<<ErrorH1<<"  "<<ErrorH1/uhH1<<"  "<<endl;
		cout<<"L2error:"<<uhL2<<"  "<<ErrorL2<<"  "<<ErrorL2/uhL2<<"  "<<endl<<endl;
//        cout<<parabolic.RHS<<endl;

		

	}
}


