#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

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

#include "vemfunction.h"

void poissonVEPkC0()
{

	for(int meshID=1;meshID<2;meshID++)
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
        
        PoissonModel poisson(VE);

		start=clock();
		poisson.Solve();
		double uhEnergy=poisson.uh.EnergyNorm(poisson.StiffMatrix);
		double EnergyError=poisson.EnergyError();
		end=clock();

		t1=double(end-start)/CLOCKS_PER_SEC;
		cout<< "meshID= "<<meshID<<endl;
		cout<<"solving time: " << t1 <<" s"<< endl;
		cout<<uhEnergy<<"  "<<EnergyError<<"  "<<EnergyError/uhEnergy<<endl<<endl;
 //       if(meshID==1) cout<<poisson.RHS<<endl;
//		for(int i=0;i<dof.Dof_Num;i++)
//			cout<<poisson.uh[i]<<"  "<<poisson.uI[i]<<endl;

	}
}


