//
//  SemilinearProbVEPkNC.h
//  Vefemz
//
//  Created by 张蓓 on 2021/3/11.
//  Copyright © 2021年 Jikun Zhao. All rights reserved.
//

#ifndef SemilinearProbVEPkNC_h
#define SemilinearProbVEPkNC_h

#include <iostream>
using namespace std;
#include "stdlib.h"
#include "stdio.h"
#include "time.h"

#include "mesh.h"

#include "virtualelement.h"
#include "problemmodel.h"

void SemilinearProbVEPkNC()
{
    FILE* fpdata=fopen("./example/SemilinearResult/error.txt","w");
 
    fprintf(fpdata,"h          ||e||         |e|_1        \n\n");
    for(int meshID=1;meshID<5;meshID++)
    {
        
        clock_t start,end;
        double t1=0;
        int k=1;
        FILE *fp;
        //        fp=GetMeshPointer(trimesh,meshID);
               fp=GetMeshPointer(rectmesh,meshID);
       // fp=GetMeshPointer(unipolymesh,meshID);
        PolyMesh ms(fp);
        fclose(fp);
        
        VEPkNC VE(k,ms);
        SemilinearEllipticModel semilinearmodel(VE);
        
        start=clock();
       semilinearmodel.Solve();
        VEMFunction error(VE),uI(VE);
        
        uI.Interpolate(semilinearmodel.pde.Solution);
        error=semilinearmodel.uh-uI;
        double uhH1norm=semilinearmodel.uh.Norm(semilinearmodel.AH1),uhL2norm=semilinearmodel.uh.Norm(semilinearmodel.AL2);
        double H1error=error.Norm(semilinearmodel.AH1),L2error=error.Norm(semilinearmodel.AL2);
        end=clock();
        
        t1=double(end-start)/CLOCKS_PER_SEC;
        cout<< "meshID= "<<meshID<<endl;
        cout<<"solving time: " << t1 <<" s"<< endl;
        cout<<"H1: "<<uhH1norm<<"  "<<H1error<<"  "<<H1error/uhH1norm<<endl;
        cout<<"L2: "<<uhL2norm<<"  "<<L2error<<"  "<<L2error/uhL2norm<<endl<<endl;
        fprintf(fpdata,"%f  %E  %E\n",ms.ElementDiameter[0],L2error,H1error);
        
    }
    
    fclose(fpdata);
}


#endif /* SemilinearProbVEPkNC_h */
