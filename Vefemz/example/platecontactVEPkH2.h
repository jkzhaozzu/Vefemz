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


void platecontactVEPkH2()
{
    FILE *fpuh;
    fpuh=fopen("./example/platecontactresult/error.txt","w");
    fprintf(fpuh,"h    EnergyError    LocalError\n");
    int k=2;//k=2
    Domain domain(-1,1,-1,1);
    FILE *fp;
   // fp=GetMeshPointer(rectmesh,7);
    fp=GetMeshPointer(unipolymesh,5);
    PolyMesh ms128(fp,domain);
    fclose(fp);
    DegreeofFreedom dof128(ms128);
    ScaledMonomialSpace sms(k);
    TriangleQuadrature TQ(k+2);
    VEPkH2 VE128(k,sms,TQ,ms128,dof128);
    PlateContactModel platecontact128(VE128,2);
    platecontact128.IterativeSolve3();
 
    
     fp=fopen("./example/platecontactresult/xyz.txt","w");
     for(int i=0;i<ms128.Node_Num;i++)
     {
     fprintf(fp,"%12.9f    %12.9f    %12.9f\n",ms128.Node[i][0],ms128.Node[i][1],platecontact128.uh[ dof128.NodeD[i][0] ]);
     }
     fclose(fp);
     
  /*
    for(int meshID=2;meshID<7;meshID++)
    {
        
        fp=GetMeshPointer(rectmesh,meshID);
        PolyMesh ms(fp,domain);
        fclose(fp);
        DegreeofFreedom dof(ms);
        VEPkH2 VE(k,sms,TQ,ms,dof);
        cout<<fixed << setprecision(5);
        PlateContactModel platecontact(VE,2);
        platecontact.IterativeSolve3();
        double uhEnergy=platecontact.uh.EnergyNorm(platecontact.StiffMatrix);
        double EnergyErr=platecontact128.EnergyError(platecontact);
        double EnergyErrFB=platecontact128.EnergyErrorFB(platecontact);
        cout<<"MeshID="<<meshID<<endl;
        cout<<"errorH2="<<EnergyErr/uhEnergy<<"     errorH2FB="<<EnergyErrFB/uhEnergy<<endl;
        
        fprintf(fpuh,"%10.8f   %10.8f    %10.8f\n",ms.ElementDiameter[0],EnergyErr/uhEnergy,EnergyErrFB/uhEnergy);
    
    }

 */
    fclose(fpuh);
}






