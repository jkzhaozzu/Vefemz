//
//  FEPkC0.cpp
//  Vefemz
//
//  Created by 张蓓 on 2019/11/10.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//
#include <iostream>
#include "finiteelement.h"

void FEPkC0::GetD(int ElemID,double **D)
{
    int i,j,EdofNum=dof.Total_Num_PerElement[ElemID];
    for(i=0;i<EdofNum;i++)//initialize
        for(j=0;j<EdofNum;j++)
            D[i][j]=0;
    
    int nv=ms.ElementVertex_Num[ElemID];
    double x[2],y[2];
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
        
    if(dof.Num_PerNode>0)
        for(i=0;i<nv;i++)
            for(j=0;j<EdofNum;j++)
                D[i][j]=SMS.GetValue(j,xpos[i],ypos[i],hE,xE,yE);
    
    GaussLobattoQuadrature GLQ(p+1);
    int *eflag=ms.ElementEdgeFlag[ElemID];
    if(dof.Num_PerEdge>0)
        for(i=0;i<nv;i++)
        {
            x[0]=xpos[i];x[1]=xpos[(i+1)%nv];
            y[0]=ypos[i];y[1]=ypos[(i+1)%nv];
            for(int k=1;k<GLQ.QuadPtsNum-1;k++)
            {
                double stdx=GLQ.QuadPts[k]; //积分点
                double jx,jy;
                if(eflag[i]==1)
                {
                    jx=(x[0]+x[1])/2.+(x[1]-x[0])*stdx/2.; jy=(y[0]+y[1])/2.+(y[1]-y[0])*stdx/2.;
                }
                else
                {
                    jx=(x[0]+x[1])/2.+(x[0]-x[1])*stdx/2.; jy=(y[0]+y[1])/2.+(y[0]-y[1])*stdx/2.;
                }
                for(int j=0;j<EdofNum;j++)
                    D[nv+i*(p-1)+(k-1)][j]=SMS.GetValue(j,jx,jy,hE,xE,yE);
            }
            
        }
    
    if(dof.Num_PerElement>0)
        for(i=0;i<dof.Num_PerElement;i++)
        {
            double vol_K=0.5*fabs((xpos[0]-xpos[2])*(ypos[1]-ypos[2])-(xpos[1]-xpos[2])*(ypos[0]-ypos[2]));
            for(int k=0;k<TQ.QuadPtsNum;k++)
            {
                double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];//积分点
                double jx=xpos[0]*j1+xpos[1]*j2+xpos[2]*j3; double jy=ypos[0]*j1+ypos[1]*j2+ypos[2]*j3;
                double jw=TQ.Weights[k];//积分权重
                for(int m=0;m<EdofNum;m++)
                    D[nv+nv*dof.Num_PerEdge+i][m]+=SMS.GetValue(m, jx, jy, hE, xE, yE)*SMS.GetValue(i, jx, jy, hE, xE, yE)*vol_K*jw/ms.ElementMeasure[ElemID];
            }
        }
    delete [] xpos; delete [] ypos;
}

void FEPkC0::GetBase(int ElemID, double **BF)
{
    int i,j,EdofNum=dof.Total_Num_PerElement[ElemID];
    double **D=new double*[EdofNum];
    VectorXd R(EdofNum),X(EdofNum);
    MatrixXd A(EdofNum,EdofNum);
    for(i=0;i<EdofNum;i++)
        D[i]=new double[EdofNum];
    GetD(ElemID,D);
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
            A(i,j)=D[i][j];
    for(i=0;i<EdofNum;i++)
    {
        for(j=0;j<EdofNum;j++)
            if(j==i)    R[j]=1;
            else    R[j]=0;
        X=A.colPivHouseholderQr().solve(R);
        for(j=0;j<EdofNum;j++)    BF[j][i]=X[j];
    }
    for(i=0;i<EdofNum;i++)  delete []D[i];
    delete[] D;
}

void FEPkC0::GetG_L2(int ElemID,double **G)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID],EdofNum=dof.Total_Num_PerElement[ElemID];
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
            G[i][j]=0;
    
    double *xpos=new double[nv];//存储单元位置
    double *ypos=new double[nv];
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    double vi,vj;
    
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
 
    double vol_K=0.5*fabs((xpos[0]-xpos[2])*(ypos[1]-ypos[2])-(xpos[1]-xpos[2])*(ypos[0]-ypos[2]));
    for(int k=0;k<TQ.QuadPtsNum;k++)
    {
        double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2]; //积分点
        double jx=xpos[0]*j1+xpos[1]*j2+xpos[2]*j3; double jy=ypos[0]*j1+ypos[1]*j2+ypos[2]*j3;
        double jw=TQ.Weights[k];//积分权重
        for(int m=0;m<EdofNum;m++)
        {
            vi=SMS.GetValue(m, jx, jy, hE, xE, yE);
            for(int n=0;n<EdofNum;n++)
            {
                vj=SMS.GetValue(n, jx, jy, hE, xE, yE);
                G[m][n]+=vi*vj*vol_K*jw;
            }
        }
    }
    delete[] xpos;  delete[] ypos;
}

void FEPkC0::GetA_L2(int ElemID,double **A)
{
    int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
    double **BF=new double*[EdofNum]; double **G=new double*[EdofNum];
    double **GB=new double*[EdofNum];
    for(i=0;i<EdofNum;i++)
    {
        BF[i]=new double[EdofNum];  G[i]=new double[EdofNum];   GB[i]=new double[EdofNum];
    }
    GetBase(ElemID, BF);  GetG_L2(ElemID, G);
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
        {
            GB[i][j]=0;
            for(k=0;k<polydim;k++)    GB[i][j]+=G[i][k]*BF[k][j];
        }
    
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
        {
            A[i][j]=0;
            for(k=0;k<EdofNum;k++)    A[i][j]+=BF[k][i]*GB[k][j];
        }
    
    for(i=0;i<EdofNum;i++)
    {
        delete []BF[i];delete []G[i];delete []GB[i];
    }
    delete []BF; delete []G; delete []GB;
}

void FEPkC0::GetG_H1(int ElemID,double **G)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID],EdofNum=dof.Total_Num_PerElement[ElemID];;
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
            G[i][j]=0;
    
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    
    double vix,viy,vjx,vjy;
    
    for(i=0;i<nv;i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];  ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    double vol_K=0.5*fabs((xpos[0]-xpos[2])*(ypos[1]-ypos[2])-(xpos[1]-xpos[2])*(ypos[0]-ypos[2]));
    for(int k=0;k<TQ.QuadPtsNum;k++)
    {
        double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];  //积分点
        double jx=xpos[0]*j1+xpos[1]*j2+xpos[2]*j3; double jy=ypos[0]*j1+ypos[1]*j2+ypos[2]*j3;
        double jw=TQ.Weights[k];  //积分权重
            
        for(int m=0;m<EdofNum;m++)
        {
            vix=SMS.GetDerivative(m,1,0, jx, jy, hE, xE, yE);
            viy=SMS.GetDerivative(m,0,1, jx, jy, hE, xE, yE);
            for(int n=0;n<EdofNum;n++)
            {
                vjx=SMS.GetDerivative(n,1,0, jx, jy, hE, xE, yE);
                vjy=SMS.GetDerivative(n,0,1, jx, jy, hE, xE, yE);
                G[m][n]+=(vix*vjx+viy*vjy)*vol_K*jw;
            }
        }
    }

    delete[] xpos;  delete[] ypos; 
}

void FEPkC0::GetA_H1(int ElemID, double **A)
{
    int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
    double **BF=new double*[EdofNum]; double **G=new double*[EdofNum];
    double **GB=new double*[EdofNum];
    for(i=0;i<EdofNum;i++)
    {
        BF[i]=new double[EdofNum];  G[i]=new double[EdofNum];   GB[i]=new double[EdofNum];
    }
    GetBase(ElemID, BF);  GetG_H1(ElemID, G);
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
        {
            GB[i][j]=0;
            for(k=0;k<polydim;k++)    GB[i][j]+=G[i][k]*BF[k][j];
        }
    
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
        {
            A[i][j]=0;
            for(k=0;k<EdofNum;k++)    A[i][j]+=BF[k][i]*GB[k][j];
        }
    
    for(i=0;i<EdofNum;i++)
    {
        delete []BF[i];delete []G[i];delete []GB[i];
    }
    delete []BF; delete []G; delete []GB;
}

void FEPkC0::GetRHS(int ElemID,FunctionP Source,double *LocF)
{
    int i,j,EdofNum=dof.Total_Num_PerElement[ElemID];
    double *LocFB=new double[EdofNum];    double **BF=new double*[EdofNum];
    for(i=0;i<EdofNum;i++)    BF[i]=new double[EdofNum];
    
    GetRHSL2B(ElemID,Source,LocFB);    GetBase(ElemID,BF);    
    for(i=0;i<EdofNum;i++)
    {
        LocF[i]=0;
        for(j=0;j<EdofNum;j++)
            LocF[i]+=LocFB[j]*BF[j][i];
    }
    
    for(i=0;i<EdofNum;i++)    delete []BF[i];
    delete []BF; delete []LocFB;
}

void FEPkC0::GetRHSL2B(int ElemID,FunctionP Source,double *LocFB)
{
    int i,j,k,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim;i++)    LocFB[i]=0;
    double *fv=new double[1];
    double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<nv;i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];   ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    double vol_K=0.5*fabs((xpos[0]-xpos[2])*(ypos[1]-ypos[2])-(xpos[1]-xpos[2])*(ypos[0]-ypos[2]));
    for(k=0;k<TQ.QuadPtsNum;k++)
    {
        double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];     //积分点
        double jx=xpos[0]*j1+xpos[1]*j2+xpos[2]*j3; double jy=ypos[0]*j1+ypos[1]*j2+ypos[2]*j3;
        double jw=TQ.Weights[k];//积分权重
        for(j=0;j<polydim;j++)
        {
            Source(jx,jy,fv);
            LocFB[j]+=fv[0]*SMS.GetValue(j,jx, jy, hE, xE, yE)*vol_K*jw;
        }
    }
    delete [] fv; delete []xpos; delete []ypos;
}

void FEPkC0::GetDofVal(FunctionP u,double *uI)
{
    int i,j,k,en,nv;
    double px,py,ex[2],ey[2];
    double value[3];
    if(dof.Num_PerNode>0)
        for(i=0;i<ms.Node_Num;i++)
        {
            px=ms.Node[i][0];    py=ms.Node[i][1];
            u(px,py,value);
            for(j=0;j<dof.Num_PerNode;j++)
                uI[ dof.NodeD[i][j] ]=value[j];
        }
    if(dof.Num_PerEdge>0)
    {
        GaussLobattoQuadrature GLQ(p+1);
        double jx,jy,stdx;
        for(i=0;i<ms.Edge_Num;i++)
        {
            en=ms.Edge[i][0];
            ex[0]=ms.Node[en][0]; ey[0]=ms.Node[en][1];
            en=ms.Edge[i][1];
            ex[1]=ms.Node[en][0]; ey[1]=ms.Node[en][1];
            
            for(j=0;j<dof.Num_PerEdge;j++)
            {
                stdx=GLQ.QuadPts[j+1]; //积分点
                jx=(ex[0]+ex[1])/2.+(ex[1]-ex[0])*stdx/2.; jy=(ey[0]+ey[1])/2.+(ey[1]-ey[0])*stdx/2.;
                u(jx,jy,value);
                uI[dof.EdgeD[i][j]]=value[0];
            }
        }
    }
    if(dof.Num_PerElement>0)
        for(int ElemID=0;ElemID<ms.Element_Num;ElemID++)
        {
            nv=ms.ElementVertex_Num[ElemID];
            double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置
            double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
            for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
            {
                xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];    ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
            }
            for(i=0;i<dof.Num_PerElement;i++)
            {
                uI[ dof.ElementD[ElemID][i] ]=0;
                double vol_K=0.5*fabs((xpos[0]-xpos[2])*(ypos[1]-ypos[2])-(xpos[1]-xpos[2])*(ypos[0]-ypos[2]));
                for(k=0;k<TQ.QuadPtsNum;k++)
                {
                    double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];//积分点
                    double jx=xpos[0]*j1+xpos[1]*j2+xpos[2]*j3; double jy=ypos[0]*j1+ypos[1]*j2+ypos[2]*j3;
                    double jw=TQ.Weights[k];//积分权重
                    u(jx,jy,value);
                    uI[ dof.ElementD[ElemID][i] ]+=value[0]*SMS.GetValue(i, jx, jy, hE, xE, yE)*vol_K*jw/ms.ElementMeasure[ElemID];
                }
            }
            delete []xpos; delete[] ypos;
        }
}

void FEPkC0::GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval)
{
    int i,j,n=0,en;
    double value[1],x,y,ex[2],ey[2];
    
    if(dof.Num_PerNode>0)
        for(i=0;i<ms.B_Node_Num;i++)
        {
            x=ms.Node[ ms.B_Node[i] ][0];    y=ms.Node[ ms.B_Node[i] ][1];
            BFunc(x,y,value);
            for(j=0;j<dof.Num_PerNode;j++)
            {
                Bdof[n]=dof.NodeD[ ms.B_Node[i] ] [j];
                Bdofval[n]=value[j];
                n++;
            }
        }
    if(dof.Num_PerEdge>0)
    {
        GaussLobattoQuadrature GLQ(p+1);
        double jx,jy,stdx;
        for(i=0;i<ms.B_Edge_Num;i++)
        {
            en=ms.Edge[ms.B_Edge[i] ][0];
            ex[0]=ms.Node[en][0]; ey[0]=ms.Node[en][1];
            en=ms.Edge[ms.B_Edge[i] ][1];
            ex[1]=ms.Node[en][0]; ey[1]=ms.Node[en][1];
            
            for(int j=0;j<dof.Num_PerEdge;j++)
            {
                Bdof[n]=dof.EdgeD[ ms.B_Edge[i] ] [j];
                stdx=GLQ.QuadPts[j+1]; //积分点
                jx=(ex[0]+ex[1])/2.+(ex[1]-ex[0])*stdx/2.; jy=(ey[0]+ey[1])/2.+(ey[1]-ey[0])*stdx/2.;
                BFunc(jx,jy,value);
                Bdofval[n]=value[0];
                n++;
            }
        }
    }
}
