//
//  mesh.cpp
//  VEM2DforMac
//
//  Created by 张蓓 on 2019/1/3.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//
#include <iostream>
#include <stdio.h>
#include "mesh.h"

PolyMesh::~PolyMesh()
{
    //节点
    int i;
    for(i=0;i<Node_Num;i++)
        delete[] Node[i];
    delete[] Node;
    delete[] NodeType;
    
    //element vertex and barycenter
    for(i=0;i<Element_Num;i++)
    {
        delete[] Element[i];
        delete[] ElementBarycenter[i];
    }
        
    delete[] Element;
    delete[] ElementMeasure;    delete[] ElementDiameter;
    delete[] ElementBarycenter;
    
    //边单元
    for(i=0;i<Element_Num;i++)
    {
        delete[] ElementEdge[i];
        delete[] ElementEdgeFlag[i];
    }
    delete[] ElementEdge;
    delete[] ElementEdgeFlag;
    
    //边
    for(i=0;i<Edge_Num;i++)
    {
        delete[] Edge[i]; delete []EdgeElement[i];
    }
    delete[] Edge;
    delete[] EdgeType;
    delete[] EdgeElement;
    
    //边界点和内部点
    delete[] B_Node;
    delete[] I_Node;
    
    //边界边和内部边
    delete[] B_Edge;
    delete[] I_Edge;
    
}

void PolyMesh::ReadMesh(FILE *fp)
{
    fscanf(fp,"%d",&Node_Num);//顶点数
    fscanf(fp,"%d",&Edge_Num);//边数
    fscanf(fp,"%d",&Element_Num);//单元数
}

void PolyMesh::MallocMeshMemory() //分配网格顶点内存空间
{
    int i;
    Node=new double*[Node_Num];//顶点
    NodeType=new int[Node_Num];
    for(i=0;i<Node_Num;i++)
        Node[i]=new double[2];
    
    Edge=new int*[Edge_Num];//边
    EdgeType=new int[Edge_Num];
    EdgeElement=new int*[Edge_Num];
    for(i=0;i<Edge_Num;i++)
    {
        Edge[i]=new int[2];EdgeElement[i]=new int[2];
        EdgeElement[i][0]=-1; EdgeElement[i][1]=-1;
    }
    
    Element=new int*[Element_Num];//单元顶点
    ElementEdge=new int*[Element_Num];//单元边
    ElementEdgeFlag=new int*[Element_Num];
    
    ElementVertex_Num=new int[Element_Num];
    ElementDiameter=new double[Element_Num];
    ElementMeasure=new double[Element_Num];
    ElementBarycenter=new double*[Element_Num];
    for(i=0;i<Element_Num;i++)
        ElementBarycenter[i]=new double[2];
}

void PolyMesh::Construct(FILE *fp)
{
    int i,m,n;
    //生成顶点
    B_Node_Num=0;
    for(i=0;i<Node_Num;i++)
    {
        fscanf(fp,"%lf%lf",&Node[i][0],&Node[i][1]);
        fscanf(fp,"%d",&NodeType[i]);
        if(NodeType[i]==1)
            B_Node_Num+=1;
    }
    I_Node_Num=Node_Num-B_Node_Num;
    
    //生成边
    B_Edge_Num=0;
    for(i=0;i<Edge_Num;i++)
    {
        fscanf(fp,"%d",&Edge[i][0]);
        fscanf(fp,"%d",&Edge[i][1]);
        fscanf(fp,"%d",&EdgeType[i]);
        if(EdgeType[i]==1)
            B_Edge_Num+=1;
    }
    I_Edge_Num=Edge_Num-B_Edge_Num;
    
    //生成单元-顶点
    for(i=0;i<Element_Num;i++)
    {
        fscanf(fp,"%d",&ElementVertex_Num[i]);
        Element[i]=new int[ElementVertex_Num[i]];
        ElementEdge[i]=new int[ElementVertex_Num[i]];
        ElementEdgeFlag[i]=new int[ElementVertex_Num[i]];
        for(int j=0;j<ElementVertex_Num[i];j++)
            fscanf(fp,"%d",&Element[i][j]);
    }
    
    //生成单元-边 和 边-单元
    for(i=0;i<Element_Num;i++)
    {
        m=ElementVertex_Num[i];
        for(int j=0;j<m;j++)
        {
            for(int k=0;k<Edge_Num;k++)
            {
                if(Element[i][j]==Edge[k][0]&&Element[i][(j+1)%m]==Edge[k][1])
                {
                    ElementEdge[i][j]=k;
                    ElementEdgeFlag[i][j]=1;
                }
                if(Element[i][j]==Edge[k][1]&&Element[i][(j+1)%m]==Edge[k][0])
                {
                    ElementEdge[i][j]=k;
                    ElementEdgeFlag[i][j]=-1;
                }
            }
        }
        
        for(int j=0;j<m;j++)
        {
            n=ElementEdge[i][j];
            if(EdgeElement[n][0]==-1)
                EdgeElement[n][0]=i;
            else
                EdgeElement[n][1]=i;
        }        
    }
    //生成边界顶点和内部顶点
    B_Node=new int[B_Node_Num];
    I_Node=new int[I_Node_Num];
    m=0; n=0;
    for(i=0;i<Node_Num;i++)
        if(NodeType[i]==1)
        {
            B_Node[m]=i; m++;
        }
        else
        {
            I_Node[n]=i; n++;
        }
    
    //生成边界边和内部边
    B_Edge=new int[B_Edge_Num];
    I_Edge=new int[I_Edge_Num];
    m=0; n=0;
    for(i=0;i<Edge_Num;i++)
        if(EdgeType[i]==1)
        {
            B_Edge[m]=i; m++;
        }
        else
        {
            I_Edge[n]=i; n++;
        }
    
    //generate element barycenter, measure and diameter
    for(i=0;i<Element_Num;i++)
    {
        double *xpos=new double[ElementVertex_Num[i]];//存储单元位置
        double *ypos=new double[ElementVertex_Num[i]];
        
        double x[3],y[3];
        int nv=ElementVertex_Num[i],j,k,m;
        for(j=0;j<nv;j++)
        {
            xpos[j]=Node[ Element[i][j] ][0];
            ypos[j]=Node[ Element[i][j] ][1];
        }
        
        ElementBarycenter[i][0]=0;  ElementBarycenter[i][1]=0;
        ElementMeasure[i]=0;    ElementDiameter[i]=0;
        for(j=0;j<nv;j++)    ElementMeasure[i]+=xpos[j]*ypos[(j+1)%nv]-xpos[(j+1)%nv]*ypos[j];
        
        for(j=0;j<nv;j++)
        {
            ElementBarycenter[i][0]+=(xpos[j]+xpos[(j+1)%nv])*(xpos[j]*ypos[(j+1)%nv]-xpos[(j+1)%nv]*ypos[j]);
            ElementBarycenter[i][1]+=(ypos[j]+ypos[(j+1)%nv])*(xpos[j]*ypos[(j+1)%nv]-xpos[(j+1)%nv]*ypos[j]);
        }
        ElementMeasure[i]=0.5*fabs(ElementMeasure[i]);
        ElementBarycenter[i][0]=ElementBarycenter[i][0]/(ElementMeasure[i]*6);
        ElementBarycenter[i][1]=ElementBarycenter[i][1]/(ElementMeasure[i]*6);
        for(j=0;j<nv;j++)
            for(k=0;k<nv-2;k++)
            {
                x[0]=xpos[j];x[1]=xpos[(j+k+1)%nv];x[2]=xpos[(j+k+2)%nv];
                y[0]=ypos[j];y[1]=ypos[(j+k+1)%nv];y[2]=ypos[(j+k+2)%nv];
                for(m=0;m<3;m++)
                {
                    double d=sqrt((x[m]-x[(m+1)%3])*(x[m]-x[(m+1)%3])+(y[(m+1)%3]-y[m])*(y[(m+1)%3]-y[m]));
                    if(ElementDiameter[i]<d)
                        ElementDiameter[i]=d;
                }
            }
    }   
}

void PolyMesh::Construct(FILE *fp,Domain &domain)
{
    int i,m,n;
    //生成顶点
    B_Node_Num=0;
    double x1=domain.x_min,x2=domain.x_max,y1=domain.y_min,y2=domain.y_max;
    double hx=x2-x1,hy=y2-y1;
    for(i=0;i<Node_Num;i++)
    {
        fscanf(fp,"%lf%lf",&Node[i][0],&Node[i][1]);
        fscanf(fp,"%d",&NodeType[i]);
        Node[i][0]=Node[i][0]*hx+x1;    Node[i][1]=Node[i][1]*hy+x1;
        
        if(NodeType[i]==1)
            B_Node_Num+=1;
    }
    I_Node_Num=Node_Num-B_Node_Num;
    
    //生成边
    B_Edge_Num=0;
    for(i=0;i<Edge_Num;i++)
    {
        fscanf(fp,"%d",&Edge[i][0]);
        fscanf(fp,"%d",&Edge[i][1]);
        fscanf(fp,"%d",&EdgeType[i]);
        if(EdgeType[i]==1)
            B_Edge_Num+=1;
    }
    I_Edge_Num=Edge_Num-B_Edge_Num;
    
    //生成单元-顶点
    for(i=0;i<Element_Num;i++)
    {
        fscanf(fp,"%d",&ElementVertex_Num[i]);
        Element[i]=new int[ElementVertex_Num[i]];
        ElementEdge[i]=new int[ElementVertex_Num[i]];
        ElementEdgeFlag[i]=new int[ElementVertex_Num[i]];
        for(int j=0;j<ElementVertex_Num[i];j++)
            fscanf(fp,"%d",&Element[i][j]);
    }
    
    //生成单元-边 和 边-单元
    for(i=0;i<Element_Num;i++)
    {
        m=ElementVertex_Num[i];
        for(int j=0;j<m;j++)
        {
            for(int k=0;k<Edge_Num;k++)
            {
                if(Element[i][j]==Edge[k][0]&&Element[i][(j+1)%m]==Edge[k][1])
                {
                    ElementEdge[i][j]=k;
                    ElementEdgeFlag[i][j]=1;
                }
                if(Element[i][j]==Edge[k][1]&&Element[i][(j+1)%m]==Edge[k][0])
                {
                    ElementEdge[i][j]=k;
                    ElementEdgeFlag[i][j]=-1;
                }
            }
        }
        for(int j=0;j<m;j++)
        {
            n=ElementEdge[i][j];
            if(EdgeElement[n][0]==-1)
                EdgeElement[n][0]=i;
            else
                EdgeElement[n][1]=i;
        }    
    }
    //生成边界顶点和内部顶点
    B_Node=new int[B_Node_Num];
    I_Node=new int[I_Node_Num];
    m=0; n=0;
    for(i=0;i<Node_Num;i++)
        if(NodeType[i]==1)
        {
            B_Node[m]=i; m++;
        }
        else
        {
            I_Node[n]=i; n++;
        }
    
    //生成边界边
    B_Edge=new int[B_Edge_Num];
    I_Edge=new int[I_Edge_Num];
    m=0; n=0;
    for(i=0;i<Edge_Num;i++)
        if(EdgeType[i]==1)
        {
            B_Edge[m]=i; m++;
        }
        else
        {
            I_Edge[n]=i; n++;
        }
    
    //generate element barycenter, measure and diameter
    for(i=0;i<Element_Num;i++)
    {
        double *xpos=new double[ElementVertex_Num[i]];//存储单元位置
        double *ypos=new double[ElementVertex_Num[i]];
        
        double x[3],y[3];
        int nv=ElementVertex_Num[i],j,k,m;
        for(j=0;j<nv;j++)
        {
            xpos[j]=Node[ Element[i][j] ][0];
            ypos[j]=Node[ Element[i][j] ][1];
        }
        
        ElementBarycenter[i][0]=0;  ElementBarycenter[i][1]=0;
        ElementMeasure[i]=0;    ElementDiameter[i]=0;
        for(j=0;j<nv;j++)    ElementMeasure[i]+=xpos[j]*ypos[(j+1)%nv]-xpos[(j+1)%nv]*ypos[j];
        
        for(j=0;j<nv;j++)
        {
            ElementBarycenter[i][0]+=(xpos[j]+xpos[(j+1)%nv])*(xpos[j]*ypos[(j+1)%nv]-xpos[(j+1)%nv]*ypos[j]);
            ElementBarycenter[i][1]+=(ypos[j]+ypos[(j+1)%nv])*(xpos[j]*ypos[(j+1)%nv]-xpos[(j+1)%nv]*ypos[j]);
        }
        ElementMeasure[i]=0.5*fabs(ElementMeasure[i]);
        ElementBarycenter[i][0]=ElementBarycenter[i][0]/(ElementMeasure[i]*6);
        ElementBarycenter[i][1]=ElementBarycenter[i][1]/(ElementMeasure[i]*6);
        for(j=0;j<nv;j++)
            for(k=0;k<nv-2;k++)
            {
                x[0]=xpos[j];x[1]=xpos[(j+k+1)%nv];x[2]=xpos[(j+k+2)%nv];
                y[0]=ypos[j];y[1]=ypos[(j+k+1)%nv];y[2]=ypos[(j+k+2)%nv];
                for(m=0;m<3;m++)
                {
                    double d=sqrt((x[m]-x[(m+1)%3])*(x[m]-x[(m+1)%3])+(y[(m+1)%3]-y[m])*(y[(m+1)%3]-y[m]));
                    if(ElementDiameter[i]<d)
                        ElementDiameter[i]=d;
                }
            }
    }
}

int PolyMesh::FindElemID(double x, double y)
{
    int i,j,nv,ElemID=-1,flag=0;
    double x0,y0,x1,y1,x2,y2;
    for (i=0;i<Element_Num;i++)
    {
        nv=ElementVertex_Num[i];
        x0=ElementBarycenter[i][0]; y0=ElementBarycenter[i][1];
        for(j=0;j<nv;j++)
        {
            x1=Node[ Element[i][j] ][0]; y1=Node[ Element[i][j] ][1];
            x2=Node[ Element[i][(j+1)%nv] ][0]; y2=Node[ Element[i][(j+1)%nv] ][1];
            if((y1-y0)*(x-x0)-(x1-x0)*(y-y0)<=0&&(y2-y1)*(x-x1)-(x2-x1)*(y-y1)<=0&&(y0-y2)*(x-x2)-(x0-x2)*(y-y2)<=0)
                ElemID=i;
        }
        if(ElemID!=-1)  break;
    }
    return ElemID;
}

FILE * GetMeshPointer(int meshtype,int meshID)
{
	FILE *fp=NULL;
	if(meshtype==trimesh)
		switch (meshID)
		{
		case 1:
			fp=fopen("./mesh/meshdata/trimesh/mesh2.dat","r");
			break;
		case 2:
			fp=fopen("./mesh/meshdata/trimesh/mesh4.dat","r");
			break;
		case 3:
			fp=fopen("./mesh/meshdata/trimesh/mesh8.dat","r");
			break;
		case 4:
			fp=fopen("./mesh/meshdata/trimesh/mesh16.dat","r");
			break;
		case 5:
			fp=fopen("./mesh/meshdata/trimesh/mesh32.dat","r");
			break;
		case 6:
			fp=fopen("./mesh/meshdata/trimesh/mesh64.dat","r");
			break;
		case 7:
			fp=fopen("./mesh/meshdata/trimesh/mesh128.dat","r");
			break;
		case 8:
			fp=fopen("./mesh/meshdata/trimesh/mesh256.dat","r");
			break;
		default:
			printf("ERROR: mesh %d不存在!\n",meshID);
			break;
		}
	else if(meshtype==rectmesh)
		switch (meshID)
		{
		case 1:
			fp=fopen("./mesh/meshdata/rectmesh/mesh2.dat","r");
			break;
		case 2:
			fp=fopen("./mesh/meshdata/rectmesh/mesh4.dat","r");
			break;
		case 3:
			fp=fopen("./mesh/meshdata/rectmesh/mesh8.dat","r");
			break;
		case 4:
			fp=fopen("./mesh/meshdata/rectmesh/mesh16.dat","r");
			break;
		case 5:
			fp=fopen("./mesh/meshdata/rectmesh/mesh32.dat","r");
			break;
		case 6:
			fp=fopen("./mesh/meshdata/rectmesh/mesh64.dat","r");
			break;
		case 7:
			fp=fopen("./mesh/meshdata/rectmesh/mesh128.dat","r");
			break;
		case 8:
			fp=fopen("./mesh/meshdata/rectmesh/mesh256.dat","r");
			break;
		default:
			printf("ERROR: mesh %d不存在!\n",meshID);
			break;
		}
	else if(meshtype==unipolymesh)
		switch (meshID)
		{
		case 1:
			fp=fopen("./mesh/meshdata/unipolymesh/mesh2.dat","r");
			break;
		case 2:
			fp=fopen("./mesh/meshdata/unipolymesh/mesh4.dat","r");
			break;
		case 3:
			fp=fopen("./mesh/meshdata/unipolymesh/mesh8.dat","r");
			break;
		case 4:
			fp=fopen("./mesh/meshdata/unipolymesh/mesh16.dat","r");
			break;
		case 5:
			fp=fopen("./mesh/meshdata/unipolymesh/mesh32.dat","r");
			break;
		case 6:
			fp=fopen("./mesh/meshdata/unipolymesh/mesh64.dat","r");
			break;
		case 7:
			fp=fopen("./mesh/meshdata/unipolymesh/mesh128.dat","r");
			break;
		default:
			printf("ERROR: mesh %d不存在!\n",meshID);
			break;
		}
	return fp;

}

