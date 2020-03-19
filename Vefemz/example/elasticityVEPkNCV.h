#include <iostream>
using namespace std;
#include "stdlib.h"
#include "stdio.h"
#include <iomanip>
#include "time.h"

#include "mesh.h"
#include "dof.h"
#include "polynomialspace.h"
#include "quadrature.h"
#include "virtualelement.h"
#include "problemmodel.h"

void elasticityVEPkNCV()
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
		ScaledMonomialSpaceV smsV(k);
		TriangleQuadrature TQ(k+2);		
		VEPkNCV VE(k,smsV,TQ,ms,dof);

		ElasticityModel elasticity(VE);
		
		start=clock();
		elasticity.Solve();
		double uhEnergy=elasticity.uh.EnergyNorm(elasticity.StiffMatrix);
		double EnergyError=elasticity.EnergyError();
		double maxerror=elasticity.MaxError();
		end=clock();

		t1=double(end-start)/CLOCKS_PER_SEC;
		cout<< "meshID= "<<meshID<<endl;
		cout<<"solving time: " << t1 <<" s"<< endl;
		cout<<uhEnergy<<" EnergyError: "<<EnergyError<<" EnergyError/uhEnergy: "<<EnergyError/uhEnergy<<" maxerror: "<<maxerror<<endl<<endl;

/*		for(int i=0;i<dof.Dof_Num;i++)
		{
			cout<<"i= "<<i<<"  "<<elasticity.uh[i]<<"  "<<elasticity.RHS[i];
			cout<<endl;
		}
*/
	}
}

