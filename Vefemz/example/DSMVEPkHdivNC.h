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
#include "mixedvirtualelement.h"
#include "problemmodel.h"



void DSMVEPkHdivNC()
{
    FILE* fpdata=fopen("./example/DSresult/error.txt","w");

    for(int meshID=2;meshID<7;meshID++)
    {
        
        clock_t start,end;
        double t1=0;
        int k=2;
        FILE *fp;
        //        fp=GetMeshPointer(trimesh,meshID);
        fp=GetMeshPointer(rectmesh,meshID);
       // fp=GetMeshPointer(unipolymesh,meshID);
        PolyMesh ms(fp);
        fclose(fp);
        DegreeofFreedom dof1(ms);
        DegreeofFreedom dof2(ms);
        ScaledMonomialSpaceV smsV(k);
        ScaledMonomialSpace sms(k-1);
        TriangleQuadrature TQ(k+2);
        VEPkHdivNCV VE1(k,smsV,TQ,ms,dof1);
        VEPkDG VE2(k-1,sms,TQ,ms,dof2);
        MVEPkHdivNC VE(VE1,VE2);

        DarcyStokesModel stokes(VE);
 
        start=clock();
        stokes.Solve();
     
//        stokes.OutputData();
        end=clock();
        double uhEnergy=stokes.uh.EnergyNorm(stokes.StiffMatrixA);
        double uhError=stokes.VelocityEnergyError();
        
        t1=double(end-start)/CLOCKS_PER_SEC;
        cout<< "meshID= "<<meshID<<endl;
        cout<<"solving time: " << t1 <<" s"<< endl;
        cout<<"uhEnergy="<<uhEnergy<<" VelocityEnergyErr=  "<<uhError<<" relative err "<<uhError/uhEnergy<<endl;
        
        double phEnergy=stokes.ph.EnergyNorm(stokes.StiffMatrixC);
        double phError=stokes.PressureL2Error();
        double maxerror=stokes.PressureMaxError();
        cout<<"phL2="<<phEnergy<<" PressureL2Err=  "<<phError<<" relative err "<<phError/uhEnergy<<endl<<endl;
        fprintf(fpdata,"%E    %E    %E\n",ms.ElementDiameter[0],uhError,phError);
 
    }
    fclose(fpdata);
}


