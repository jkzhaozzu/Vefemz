//
//  mesh.h
//  VEM2DforMac
//
//  Created by Jikun Zhao on 2019/1/2.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//

#ifndef mesh_h
#define mesh_h

#include "assert.h"
#include "math.h"
#include "stdio.h"

//domain class
class Domain
{
public:
    Domain(){};
    Domain(double x1,double x2,double y1,double y2)
    {
        x_min=x1; x_max=x2; y_min=y1; y_max=y2;
    };
    ~Domain(){};
    double x_min,x_max,y_min,y_max;
};

//mesh class
class PolyMesh
{

public:
    PolyMesh(){;};
	PolyMesh(FILE *fp)
	{
		ReadMesh(fp);
		MallocMeshMemory();
		Construct(fp);
	};
    PolyMesh(FILE *fp,Domain &dom):domain(dom)
    {
        ReadMesh(fp);
        MallocMeshMemory();
        Construct(fp,domain);
    };
    ~PolyMesh();
    void ReadMesh(FILE *fp);//read mesh data from file
    void MallocMeshMemory();//分配网格顶点内存空间
    void Construct(FILE *fp);//construct mesh structure
    void Construct(FILE *fp,Domain &domain);//construct mesh structure on domain (a1,b1)*(a2,b2)
    int FindElemID(double x,double y); //find mesh element ID in which the point (x,y) is
       
    //mesh information
    Domain domain;
    int Node_Num;//节点数
    int Edge_Num;//边数
	int Element_Num;//单元数
    int *ElementVertex_Num;//单元顶点数
    int B_Node_Num;//边界点数
    int B_Edge_Num;//边界边数
//    int B_Element_Num;//边界单元数
    int I_Node_Num;//内部点数
    int I_Edge_Num;//内部边数
    //mesh structure
    double **Node;//节点地址
    int *NodeType;//节点类型, 1表示边界点，0表示内部点
	int **Edge;//边地址
    int *EdgeType;//边类型, 1表示边界边，0表示内部边
    int **Element;//单元节点地址, 1表示边界单元，0表示内部单元
//    int*ElementType;//单元类型，
    int **ElementEdge;//单元边地址
    int **EdgeElement;//边单元地址，-1表示空地址
    int **ElementEdgeFlag;//单元边标志=1,正向；=-1 反向
    int *B_Node;//边界点地址
    int *B_Edge;//边界边地址
    int *I_Node;//内部点地址
    int *I_Edge;//内部边地址
//    int *B_Element;//边界单元地址
       
    double *ElementDiameter;
    double *ElementMeasure;
    double ** ElementBarycenter;    
};



enum{trimesh, rectmesh, unipolymesh};
FILE * GetMeshPointer(int meshtype,int meshID);




#endif /* mesh_h */
