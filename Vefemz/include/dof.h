//
//  dofclass.h
//  VEM2DforMac
//
//  Created by 张蓓 on 2019/1/2.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//

#ifndef dofclass_h
#define dofclass_h

#include "mesh.h"

//自由度
typedef struct Dof
{
    double value;
    int type;//记录自由度类型
    int id;//自由度所在对象号
}Dof;

//自由度类型
enum{Node_Dof,Edge_Dof,Element_Dof};

class DegreeofFreedom
{
public:
    
    DegreeofFreedom();
    DegreeofFreedom(PolyMesh &mesh):ms(mesh){;};
    ~DegreeofFreedom();
    void GetDofInfo();//根据网格信息和自由度信息计算总自由度数
    void MallocDofMemory();//分配自由度内存空间
    void Construct();
    
    //mesh information
    PolyMesh &ms;
    
    // DoF information
    int Num_PerNode;//每个顶点自由度数
    int Num_PerElement;//number of each element's interior dof 
    int Num_PerEdge;//每条边自由度数
    int *Total_Num_PerElement;//number of each element's total dof
    int Dof_Num;//自由度总个数
    
    // DoF structure
    Dof* DList;//自由度
    short *DType;//0 for inner dof and 1 for boundary dof
    int **ElementD;//each element's interior dof ID
    int **EdgeD;//每条边自由度号
    int **NodeD;//每个顶点自由度号
    int **TotalD;//each element's total dof IDs
    
};


#endif /* dofclass_h */
