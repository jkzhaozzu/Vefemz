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



void StokesMVEPkHdivNC()
{
    
    for(int meshID=1;meshID<7;meshID++)
    {
        
        clock_t start,end;
        double t1=0,t2=0;
        int k=1;
        FILE *fp;
        //        fp=GetMeshPointer(trimesh,meshID);
                fp=GetMeshPointer(rectmesh,meshID);
        //fp=GetMeshPointer(unipolymesh,meshID);
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
    
        StokesModel stokes(VE);
        
        start=clock();
        stokes.Solve();
        end=clock();
        double uhEnergy=stokes.uh.EnergyNorm(stokes.StiffMatrixA);
        
        double EnergyError=stokes.VelocityEnergyError();
        
        t1=double(end-start)/CLOCKS_PER_SEC;
        cout<< "meshID= "<<meshID<<endl;
        cout<<"solving time: " << t1 <<" s"<< endl;
        cout<<"uhEnergy="<<uhEnergy<<" VelocityEnergyErr=  "<<EnergyError<<" relative err "<<EnergyError/uhEnergy<<endl;
        
        uhEnergy=stokes.ph.EnergyNorm(stokes.StiffMatrixC);
        EnergyError=stokes.PressureL2Error();
        double maxerror=stokes.PressureMaxError();
        cout<<"phL2="<<uhEnergy<<" PressureL2Err=  "<<EnergyError<<" relative err "<<EnergyError/uhEnergy<<endl<<endl;
        
        
        //        cout<<"pI.aver="<<stokes.pI.GetAverageValue()<<" ph.aver="<<stokes.ph.GetAverageValue()<<endl;
        //        for(int i=0;i<stokes.dof2.Dof_Num;i++)
        //            cout<<"i= "<<i<<"  "<<stokes.pI[i]<<"  "<<stokes.ph[i]<<"  "<<endl;
    }
}



