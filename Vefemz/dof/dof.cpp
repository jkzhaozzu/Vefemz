//
//  dof.cpp
//  VEM2DforMac
//
//  Created by 张蓓 on 2019/1/3.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//

#include <stdio.h>
#include "dof.h"
#include "mesh.h"

DegreeofFreedom::~DegreeofFreedom()
{
    delete[] DList;
    if(Num_PerElement>0)
    {
        for(int i=0;i<ms.Element_Num;i++)
        {
            delete[] ElementD[i];
        }
        
        delete[] ElementD;
    }
    
    if(Num_PerEdge>0)
    {
        for(int i=0;i<ms.Edge_Num;i++)
            delete[] EdgeD[i];
        delete[] EdgeD;
    }
    
    if(Num_PerNode>0)
    {
        for(int i=0;i<ms.Node_Num;i++)
            delete[] NodeD[i];
        delete[] NodeD;
    }
    
    for(int i=0;i<ms.Element_Num;i++)
        delete[] TotalD[i];
    delete[] TotalD;
    delete[] DType;
}

void DegreeofFreedom::GetDofInfo()
{

    Dof_Num=Num_PerEdge*ms.Edge_Num+Num_PerElement*ms.Element_Num+Num_PerNode*ms.Node_Num;
    Total_Num_PerElement=new int[ms.Element_Num];
    for(int i=0;i<ms.Element_Num;i++)
        Total_Num_PerElement[i]=ms.ElementVertex_Num[i]*Num_PerNode+ms.ElementVertex_Num[i]*Num_PerEdge+Num_PerElement;
}

void DegreeofFreedom::MallocDofMemory()
{
    DList=new Dof[Dof_Num];
    DType=new short[Dof_Num];
    
    if(Num_PerNode>0)
    {
        NodeD=new int*[ms.Node_Num];
        for(int i=0;i<ms.Node_Num;i++)
            NodeD[i]=new int[Num_PerNode];
    }
    
    if(Num_PerEdge>0)
    {
        EdgeD=new int*[ms.Edge_Num];
        for(int i=0;i<ms.Edge_Num;i++)
            EdgeD[i]=new int[Num_PerEdge];
    }
    
    if(Num_PerElement>0)
    {
        ElementD=new int*[ms.Element_Num];
        for(int i=0;i<ms.Element_Num;i++)
            ElementD[i]=new int[Num_PerElement];
    }
    TotalD=new int*[ms.Element_Num];
    for(int i=0;i<ms.Element_Num;i++)
        TotalD[i]=new int[Total_Num_PerElement[i]];
}

void DegreeofFreedom::Construct()
{
    //产生自由度，由于要使用压缩存储，故不必考虑带宽
    int i,j,n=0;
    int tn=0;//记录总自由度数
    if(Num_PerNode>0)
    {
        for(i=0;i<ms.Node_Num;i++)
        {
            for(j=0;j<Num_PerNode;j++)
            {
                NodeD[i][j]=tn;//存储节点自由度号
                //存储自由度类型
                DList[tn].type=Node_Dof;
                DType[tn]=0;
                DList[tn++].id=i;
            }
        }
    }
    
    if(Num_PerEdge>0)
    {
        for(i=0;i<ms.Edge_Num;i++)
        {
            for(j=0;j<Num_PerEdge;j++)
            {
                EdgeD[i][j]=tn;//存储节点自由度号
                
                //存储自由度类型
                DList[tn].type=Edge_Dof;
                DType[tn]=0;
                DList[tn++].id=i;
            }
        }
    }
    
    if(Num_PerElement>0)
    {
        for(i=0;i<ms.Element_Num;i++)
        {
            for(j=0;j<Num_PerElement;j++)
            {
                ElementD[i][j]=tn;//存储节点自由度号
                
                //存储自由度类型
                DList[tn].type=Element_Dof;
                DType[tn]=0;
                DList[tn++].id=i;
            }
        }
    }
    
    //装配每个有限元单元的自由度
    for(i=0;i<ms.Element_Num;i++)
    {
        int pn=Total_Num_PerElement[i];
        int dn=0;
        for(j=0;j<ms.ElementVertex_Num[i];j++)
        {
            int nid=ms.Element[i][j];
            for(int k=0;k<Num_PerNode;k++)
                TotalD[i][dn++]=NodeD[nid][k];
        }
        for(j=0;j<ms.ElementVertex_Num[i];j++)
        {
//            int nid=ms.Element[i][j];
            int eid=ms.ElementEdge[i][j];
            for(int k=0;k<Num_PerEdge;k++)
                TotalD[i][dn++]=EdgeD[eid][k];
        }
        
        for(int k=0;k<Num_PerElement;k++)
            TotalD[i][dn++]=ElementD[i][k];
        
        assert(dn==pn);
    }
    
    if(Num_PerNode>0)
        for(i=0;i<ms.B_Node_Num;i++)
            for(j=0;j<Num_PerNode;j++)
            {
                n=NodeD[ ms.B_Node[i] ] [j];
                DType[n]=1;
            }
    
    if(Num_PerEdge>0)
        for(i=0;i<ms.B_Edge_Num;i++)
            for(j=0;j<Num_PerEdge;j++)
            {
                n=EdgeD[ ms.B_Edge[i] ] [j];
                DType[n]=1;
            }
}
