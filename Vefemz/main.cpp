//
//  main.cpp
//  VEM2DforMac
//
//  Created by Jikun Zhao on 2019/1/2.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//

#include <iostream>
#include "stdlib.h"
#include "stdio.h"
using namespace std;

//#include "./example/poissonVEPkC0.h"
//#include "./example/poissonVEPkNC.h"
//#include "./example/ParabolicVEPkC0.h"
//#include "./example/ParabolicVEPkNC.h"
//#include "./example/elasticityVEPkNCV.h"
//#include "./example/StokesMVEPkNC.h"
//#include "./example/StokesMVEPkHdivNC.h"
//#include "./example/NavierStokesMVEPkNC.h"
//#include "./example/BiharmVEPkH2.h"
//#include "./example/BiharmVEC0PkH2NC.h"
//#include "./example/BiharmVEMorleyType.h"
//#include "./example/platecontactVEC0PkH2NC.h"
//#include "./example/platecontMorleyType.h"
//#include "./example/platecontactVEPkH2.h"
#include "./example/DSMVEPkHdivNC.h"

//#include "./FEexample/poissonFEPkC0.h"
//#include "./FEexample/elasticityMFEP1NCRS.h"



int main(int argc, const char * argv[])
{

    
    DSMVEPkHdivNC();
  //  std::vector<Point> l_points = {Point(10,0), Point(20,0), Point(20,20), Point(0,20), Point(0,10), Point(10,10)};

    cout << "Hello World!\n" <<endl;
    getchar();
    return 0;
    
}

