//
//  CookMembSVEPkNC.h
//  Vefemz
//
//  Created by 张蓓 on 12/24/21.
//  Copyright © 2021 Jikun Zhao. All rights reserved.
//

#ifndef CookMembSVEPkNC_h
#define CookMembSVEPkNC_h

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

void CookMembSVEPkNC()
{
    for(int meshID=1;meshID<8;meshID++)
    {
        clock_t start,end;
        int k=1;
        FILE *fp;
        //        fp=GetMeshPointer(trimesh,meshID);
        fp=GetMeshPointer(rectmesh,meshID);
        //        fp=GetMeshPointer(unipolymesh,meshID);
        
        PolyMesh ms(fp);
        fclose(fp);
        VEPkNCV VE(k,ms);
        
        CookMembraneModel elasticity(VE);
        
        start=clock();
        elasticity.Solve();
        double value[2];
        
    //    cout<<value[0]<<" "<<value[1]<<endl;
    //    cout<<elasticity.RHS<<endl;
    //    cout<<elasticity.pde.mu<<" "<<elasticity.pde.lambda<<endl;
        VEMFunction uI(VE);
       
        double uhEnergy=elasticity.uh.Norm(elasticity.AH1);
        cout<< "meshID= "<<meshID<<endl;
        elasticity.uh.Projection(48,60,value,"H1");
        cout<<value[0]<<" "<<value[1]<<endl;

        
  
       


    }
    

}

#endif /* CookMembSVEPkNC_h */
