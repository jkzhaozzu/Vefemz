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
#include "finiteelement.h"
#include "problemmodel.h"
#include "FEmodel.h"

void poissonFEPkC0()
{
    
    for(int meshID=1;meshID<7;meshID++)
    {
        
        clock_t start,end;
        double t1=0;
        int k=1;
        FILE *fp;
        fp=GetMeshPointer(trimesh,meshID);
        PolyMesh ms(fp);
        fclose(fp);
        DegreeofFreedom dof(ms);
        ScaledMonomialSpace sms(k);
        TriangleQuadrature TQ(k+2);
        FEPkC0 FE(k,sms,TQ,ms,dof);

        FEPoissonModel poisson(FE);

        start=clock();
        poisson.Solve();
 //       cout<<poisson.StiffMatrix<<endl;
        double uhEnergy=poisson.uh.EnergyNorm(poisson.StiffMatrix);
        double EnergyError=poisson.EnergyError();
        end=clock();
        
        t1=double(end-start)/CLOCKS_PER_SEC;
        cout<< "meshID= "<<meshID<<endl;
        cout<<"solving time: " << t1 <<" s"<< endl;
        cout<<uhEnergy<<"  "<<EnergyError<<"  "<<EnergyError/uhEnergy<<endl<<endl;
    }
}



