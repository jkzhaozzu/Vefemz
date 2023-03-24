//
//  ElasticitySVEPkNCV.h
//  Vefemz
//
//  Created by 张蓓 on 2021/4/19.
//  Copyright © 2021年 Jikun Zhao. All rights reserved.
//

#ifndef ElasticitySVEPkNCV_h
#define ElasticitySVEPkNCV_h

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

void ElasticitySVEPkNCV()
{
    for(int meshID=2;meshID<5;meshID++)
    {
        clock_t start,end;
        double t1=0;
        int k=1;
        FILE *fp;
        //        fp=GetMeshPointer(trimesh,meshID);
        fp=GetMeshPointer(rectmesh,meshID);
       //         fp=GetMeshPointer(unipolymesh,meshID);
        
        PolyMesh ms(fp);
        fclose(fp);
        VEPkNCV VE(k,ms);
        
        ElasticityStabilizationModel elasticity(VE);
        
        start=clock();
        elasticity.Solve();
        
        VEMFunction error(VE),uI(VE);
        uI.Interpolate(elasticity.pde.Solution);
        error=elasticity.uh-uI;
        double uhH1=elasticity.uh.Norm(elasticity.AH1),uhL2=elasticity.uh.Norm(elasticity.AL2);
        double H1Error=error.Norm(elasticity.AH1),L2Error=error.Norm(elasticity.AL2);
        end=clock();
        
        t1=double(end-start)/CLOCKS_PER_SEC;
        cout<< "meshID= "<<meshID<<endl;
        cout<<"solving time: " << t1 <<" s"<< endl;
        cout<<"uhH1:"<<uhH1<<" uhH1Error: "<<H1Error<<" uhH1Error/uhH1: "<<H1Error/uhH1<<endl;
        cout<<"uhL2:"<<uhL2<<" uhL2Error: "<<L2Error<<" uhL2Error/uhL2: "<<L2Error/uhL2<<endl<<endl;

        //    cout << elasticity.uh << endl;
     //   cout << elasticity.StiffMatrix<<endl;
    }
    

}




#endif /* ElasticitySVEPkNCV_h */
