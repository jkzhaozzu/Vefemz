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


void BiharmVEPkH2()
{
    for(int meshID=1;meshID<7;meshID++)
    {
        
        clock_t start,end;
        double t1=0;
        int k=4;//k>=2
        FILE *fp;
      //  fp=GetMeshPointer(trimesh,meshID);
        fp=GetMeshPointer(rectmesh,meshID);
      //  fp=GetMeshPointer(unipolymesh,meshID);
        PolyMesh ms(fp);
        fclose(fp);
        DegreeofFreedom dof(ms);
        ScaledMonomialSpace sms(k);
        TriangleQuadrature TQ(k+2);
        VEPkH2 VE(k,sms,TQ,ms,dof);
        
        cout<<fixed << setprecision(5);

        BiharmonicModel biharmonicequation(VE);
        
        start=clock();
        biharmonicequation.Solve2();
        double uhEnergy=biharmonicequation.uh.EnergyNorm(biharmonicequation.StiffMatrix);
        double EnergyError=biharmonicequation.EnergyError();
        double MaxError=biharmonicequation.MaxError();
        end=clock();
        
        t1=double(end-start)/CLOCKS_PER_SEC;
        cout<< "meshID= "<<meshID<<endl;
        cout<<"solving time: " << t1 <<" s"<< endl;
        cout<<uhEnergy<<"  "<<MaxError<<"  "<<EnergyError/uhEnergy<<endl<<endl;
      //  cout<<biharmonicequation.RHS<<endl;
        //        for(int i=0;i<dof.Dof_Num;i++)
        //            cout<<"i="<<i<<" "<<biharmonicequation.uh[i]<<"  "<<biharmonicequation.RHS[i]<<endl;
        
    }
}



