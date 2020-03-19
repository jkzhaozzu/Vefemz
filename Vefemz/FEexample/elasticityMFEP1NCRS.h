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
#include "mixedfiniteelement.h"
#include "FEmodel.h"

void elasticityMFEP1NCRS()
{    
    for(int meshID=1;meshID<7;meshID++)
    {
        
        clock_t start,end;
        double t1=0;
        int k=1;
        Domain domain(-1,1,-1,1);
        FILE *fp;
        fp=GetMeshPointer(rectmesh,meshID);
        PolyMesh ms(fp,domain);
        fclose(fp);
        DegreeofFreedom dof1(ms),dof2(ms);
        TriangleQuadrature TQ(6);
        MFEP1NCRS MFE(TQ,ms,dof1,dof2);
        
        FEElasticityMixedModel elasticity(MFE);
        start=clock();
        elasticity.Solve();

        double uhH1norm=elasticity.uhH1Norm(),uhMaxerror=elasticity.uhMaxError();
        double uhL2norm=elasticity.uhL2Norm();
        double sigmahL2norm=elasticity.sigmahL2Norm(),sigmahL2error=elasticity.sigmahL2Error(),
        sigmahMaxeror=elasticity.sigmahMaxError();
        double uhH1error=elasticity.uhH1Error(),uhL2error=elasticity.uhL2Error();
        end=clock();
        
        t1=double(end-start)/CLOCKS_PER_SEC;
        cout<< "meshID= "<<meshID<<endl;
        cout<<"solving time: " << t1 <<" s "<<endl;
        cout<<"sigma L2: "<<sigmahL2norm<<"  "<<sigmahL2error<<"  "<<sigmahL2error/sigmahL2norm<<endl;
        cout<<"u H1: "<<uhH1norm<<"  "<<uhH1error<<"  "<<uhH1error/uhH1norm<<endl;
        cout<<"u L2: "<<uhL2norm<<"  "<<uhL2error<<"  "<<uhL2error/uhL2norm<<endl<<endl;



    }

}




