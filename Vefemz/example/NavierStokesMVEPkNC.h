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
#include "mixedvirtualelement.h"
#include "problemmodel.h"



void NavierStokesMVEPkNC()
{
//	FILE *fp1=fopen("error.txt","w");
//	fprintf(fp1,"Mesh#   |u_I-u_h|_1  ||p-p_h||\n");

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
		DegreeofFreedom dof1(ms);
		DegreeofFreedom dof2(ms);
		ScaledMonomialSpaceV smsV(k);
		ScaledMonomialSpace sms(k-1);
		TriangleQuadrature TQ(k+2);
        VEPkNCV VE1(k,smsV,TQ,ms,dof1);
        VEPkDG VE2(k-1,sms,TQ,ms,dof2);
        MVEPkNC MVE(VE1,VE2);

		NavierStokesModel navierstokes(MVE);
		
		start=clock();
		navierstokes.Solve();
		end=clock();
		double uhEnergy=navierstokes.uh.EnergyNorm(navierstokes.MatrixA);
		double uhEnergyError=navierstokes.VelocityEnergyError();

		t1=double(end-start)/CLOCKS_PER_SEC;
		cout<< "meshID= "<<meshID<<endl;
		cout<<"solving time: " << t1 <<" s"<< endl;
		cout<<"uhEnergy="<<uhEnergy<<" VelocityEnergyErr=  "<<uhEnergyError<<" relative err "<<uhEnergyError/uhEnergy<<endl;

		double phEnergy=navierstokes.ph.EnergyNorm(navierstokes.MatrixC);
		double phEnergyError=navierstokes.PressureL2Error();
		cout<<"phL2="<<phEnergy<<" PressureL2Err=  "<<phEnergyError<<" relative err "<<phEnergyError/phEnergy<<endl<<endl;

//		fprintf(fp1,"%d   %f   %f\n",meshID,uhEnergyError,phEnergyError);
//        cout<<navierstokes.StiffMatrix<<endl;
 //       cout<<navierstokes.MatrixAdivT<<endl;

//		cout<<"pI.aver="<<stokes.pI.GetAverageValue()<<" ph.aver="<<stokes.ph.GetAverageValue()<<endl;
//		for(int i=0;i<stokes.dof2.Dof_Num;i++)
//			cout<<"i= "<<i<<"  "<<stokes.pI[i]<<"  "<<stokes.ph[i]<<"  "<<endl;
	}
//	fclose(fp1);
}

