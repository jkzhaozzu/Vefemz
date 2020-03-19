//
//  MFEP1NCRS.cpp
//  Vefemz
//
//  Created by 张蓓 on 2020/2/17.
//  Copyright © 2020年 Jikun Zhao. All rights reserved.
//

#include <iostream>
#include "mixedfiniteelement.h"

void MFEP1NCRS::GetBasisFuntion1(int ElemID,int i,double x,double y,double *value)
{
    double xE,yE,hE;
    xE=ms.ElementBarycenter[ElemID][0];yE=ms.ElementBarycenter[ElemID][1];hE=ms.ElementDiameter[ElemID];
    switch(i)
    {
        case 0:
            value[0]=1;value[1]=0;value[2]=0; //sigma11,sigma12,sigma22
            value[3]=0;value[4]=0;//div(sigma)
            break;
        case 1:
            value[0]=3*(x-xE)*2*sqrt(2)/hE;value[1]=0;value[2]=0;
            value[3]=3*2*sqrt(2)/hE;value[4]=0;
            break;
        case 2:
            value[0]=0;value[1]=1;value[2]=0;
            value[3]=0;value[4]=0;
            break;
        case 3:
            value[0]=0;value[1]=0;value[2]=1;
            value[3]=0;value[4]=0;
            break;
        case 4:
            value[0]=0;value[1]=0;value[2]=3*(y-yE)*2*sqrt(2)/hE;
            value[3]=0;value[4]=3*2*sqrt(2)/hE;
            break;
        default:
            printf("ERROR:形函数%d不存在!\n",i);
    }
}

void MFEP1NCRS::GetBasisFuntion2(int ElemID,int i,double x,double y,double *value)
{
    double xE,yE,hE;
    xE=ms.ElementBarycenter[ElemID][0];yE=ms.ElementBarycenter[ElemID][1];hE=ms.ElementDiameter[ElemID];
    double aa=2*sqrt(2)/hE;
    switch(i)
    {
        case 0:
            value[0]=0.75-0.5*(y-yE)*aa-0.75*(x-xE)*(x-xE)*aa*aa;value[1]=0;//(u1,u2)
            value[2]=-1.5*(x-xE)*aa*aa;value[3]=-0.5*aa; //grad(u1)
            value[4]=0;value[5]=0;//grad(u2)
            break;
        case 1:
            value[0]=0;value[1]=-0.25-0.5*(y-yE)*aa+0.75*(y-yE)*(y-yE)*aa*aa;
            value[2]=0;value[3]=0;
            value[4]=0;value[5]=-0.5*aa+1.5*(y-yE)*aa*aa;
            break;
        case 2:
            value[0]=-0.25+0.5*(x-xE)*aa+0.75*(x-xE)*(x-xE)*aa*aa;value[1]=0;
            value[2]=0.5*aa+1.5*(x-xE)*aa*aa;value[3]=0;
            value[4]=0;value[5]=0;
            break;
        case 3:
            value[0]=0;value[1]=0.75+0.5*(x-xE)*aa-0.75*(y-yE)*(y-yE)*aa*aa;
            value[2]=0;value[3]=0;
            value[4]=0.5*aa;value[5]=-1.5*(y-yE)*aa*aa;
            break;
        case 4:
            value[0]=0.75+0.5*(y-yE)*aa-0.75*(x-xE)*(x-xE)*aa*aa;value[1]=0;
            value[2]=-1.5*(x-xE)*aa*aa;value[3]=0.5*aa;
            value[4]=0;value[5]=0;
            break;
        case 5:
            value[0]=0;value[1]=-0.25+0.5*(y-yE)*aa+0.75*(y-yE)*(y-yE)*aa*aa;
            value[2]=0;value[3]=0;
            value[4]=0;value[5]=0.5*aa+1.5*(y-yE)*aa*aa;
            break;
        case 6:
            value[0]=-0.25-0.5*(x-xE)*aa+0.75*(x-xE)*(x-xE)*aa*aa;value[1]=0;
            value[2]=-0.5*aa+1.5*(x-xE)*aa*aa;value[3]=0;
            value[4]=0;value[5]=0;
            break;
        case 7:
            value[0]=0;value[1]=0.75-0.5*(x-xE)*aa-0.75*(y-yE)*(y-yE)*aa*aa;
            value[2]=0;value[3]=0;
            value[4]=-0.5*aa;value[5]=-1.5*(y-yE)*aa*aa;
            break;
        default:
            printf("ERROR:形函数%d不存在!\n",i);
    }
}

void MFEP1NCRS::GetMixedAL2(int ElemID,double **A)
{
    int i,j,k,m,n;
    for(i=0;i<polydim1;i++)
        for(j=0;j<polydim1;j++)
            A[i][j]=0;
    int nv=ms.ElementVertex_Num[ElemID];
    double x[3],y[3],v1[5],v2[5];
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    double xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    for(j=0;j<nv;j++)
    {
        x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
        y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];//积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            for(m=0;m<polydim1;m++)
            {
                GetBasisFuntion1(ElemID, m, jx, jy, v1);
                for(n=0;n<polydim1;n++)
                {
                    GetBasisFuntion1(ElemID, n, jx, jy, v2);
                    A[m][n]+=(v1[0]*v2[0]+2*v1[1]*v2[1]+v1[2]*v2[2])*vol_K*jw;
                }
            }
        }
    }
//    A[0][1]=0;A[1][0]=0;A[3][4]=0;A[4][3]=0;
    delete [] xpos; delete [] ypos;
}

void MFEP1NCRS::GetMixedAL2trace(int ElemID,double **A)
{
    int i,j,k,m,n;
    for(i=0;i<polydim1;i++)
        for(j=0;j<polydim1;j++)
            A[i][j]=0;
    int nv=ms.ElementVertex_Num[ElemID];
    double x[3],y[3],v1[5],v2[5];
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    double xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    for(j=0;j<nv;j++)
    {
        x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
        y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];//积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            for(m=0;m<polydim1;m++)
            {
                GetBasisFuntion1(ElemID, m, jx, jy, v1);
                for(n=0;n<polydim1;n++)
                {
                    GetBasisFuntion1(ElemID, n, jx, jy, v2);
                    A[m][n]+=(v1[0]+v1[2])*(v2[0]+v2[2])*vol_K*jw;
                }
            }
        }
    }
    delete [] xpos; delete [] ypos;
}

void MFEP1NCRS::GetMixedAdiv(int ElemID,double **A)
{
    int i,j;
    double v1[5],v2[5],vol=ms.ElementMeasure[ElemID],hE=ms.ElementDiameter[ElemID];
    for(i=0;i<polydim1;i++)
    {
        GetBasisFuntion1(ElemID, i, 0, 0, v1);
        for(j=0;j<polydim1;j++)
        {
            GetBasisFuntion1(ElemID, j, 0, 0, v2);
            A[i][j]=hE*hE*(v1[3]*v2[3]+v1[4]*v2[4])*vol;
        }
    }
}

void MFEP1NCRS::GetMixedB(int ElemID,double **B)
{
    int i,j,k,m,n;
    for(i=0;i<polydim1;i++)
        for(j=0;j<polydim2;j++)
            B[i][j]=0;
    int nv=ms.ElementVertex_Num[ElemID];
    double x[3],y[3],v1[5],v2[6];
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    double xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    for(j=0;j<nv;j++)
    {
        x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
        y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];//积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            for(m=0;m<polydim1;m++)
            {
                GetBasisFuntion1(ElemID, m, jx, jy, v1);
                for(n=0;n<polydim2;n++)
                {
                    GetBasisFuntion2(ElemID, n, jx, jy, v2);
                    B[m][n]-=(v1[0]*v2[2]+v1[1]*(v2[3]+v2[4])+v1[2]*v2[5])*vol_K*jw;
                }
            }
        }
    }
    delete [] xpos; delete [] ypos;
}

void MFEP1NCRS::GetMixedC1(int ElemID,double **C)
{
    int i,j,k,m;
    for(i=0;i<polydim2;i++)
        for(j=0;j<polydim2;j++)
            C[i][j]=0;
    int nv=ms.ElementVertex_Num[ElemID];
    double x[2],y[2],v1[6],v2[6];
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    GaussLobattoQuadrature GLQ(4);
    for(i=0;i<nv;i++)
    {
        x[0]=xpos[i];x[1]=xpos[(i+1)%nv];
        y[0]=ypos[i];y[1]=ypos[(i+1)%nv];
        for(k=0;k<GLQ.QuadPtsNum;k++)
        {
            double stdx=GLQ.QuadPts[k]; //积分点
            double jx,jy;
            jx=(x[0]+x[1])/2.+(x[1]-x[0])*stdx/2.; jy=(y[0]+y[1])/2.+(y[1]-y[0])*stdx/2.;
            for(j=0;j<polydim2;j++)
            {
                GetBasisFuntion2(ElemID, j, jx, jy, v1);
                for(m=0;m<polydim2;m++)
                {
                    GetBasisFuntion2(ElemID, m, jx, jy, v2);
                    C[j][m]+=(v1[0]*v2[0]+v1[1]*v2[1])*GLQ.Weights[k]/2;
                }
            }
        }
    }
    delete [] xpos; delete [] ypos;    
}

void MFEP1NCRS::GetMixedC2(int EdgeID,double **C)
{
    int i,j,k;
    for(i=0;i<polydim2;i++)
        for(j=0;j<polydim2;j++)
            C[i][j]=0;
    int E1ID=ms.EdgeElement[EdgeID][0],E2ID=ms.EdgeElement[EdgeID][1];
    double x[2],y[2],v1[6],v2[6];
    x[0]=ms.Node[ ms.Edge[EdgeID][0] ][0];x[1]=ms.Node[ ms.Edge[EdgeID][1] ][0];
    y[0]=ms.Node[ ms.Edge[EdgeID][0] ][1];y[1]=ms.Node[ ms.Edge[EdgeID][1] ][1];
    GaussLobattoQuadrature GLQ(3);
    if(E2ID!=-1)
    {
        for(k=0;k<GLQ.QuadPtsNum;k++)
        {
            double stdx=GLQ.QuadPts[k]; //积分点
            double jx,jy;
            jx=(x[0]+x[1])/2.+(x[1]-x[0])*stdx/2.; jy=(y[0]+y[1])/2.+(y[1]-y[0])*stdx/2.;
            for(i=0;i<polydim2;i++)
            {
                GetBasisFuntion2(E1ID, i, jx, jy, v1);
                for(j=0;j<polydim2;j++)
                {
                    GetBasisFuntion2(E2ID, j, jx, jy, v2);
                    C[i][j]-=(v1[0]*v2[0]+v1[1]*v2[1])*GLQ.Weights[k]/2;
                }
            }
        }
    }
}

/*void MFEP1NCRS::GetMixedC3(int EdgeID,double **C)
{
    int i,j,k;
    for(i=0;i<polydim2;i++)
        for(j=0;j<polydim2;j++)
            C[i][j]=0;
    int E1ID=ms.EdgeElement[EdgeID][0],E2ID=ms.EdgeElement[EdgeID][1];
    double x[2],y[2],v1[6],v2[6];
    x[0]=ms.Node[ ms.Edge[EdgeID][0] ][0];x[1]=ms.Node[ ms.Edge[EdgeID][1] ][0];
    y[0]=ms.Node[ ms.Edge[EdgeID][0] ][1];y[1]=ms.Node[ ms.Edge[EdgeID][1] ][1];
    GaussLobattoQuadrature GLQ(3);
    if(E2ID!=-1)
    {
        for(k=0;k<GLQ.QuadPtsNum;k++)
        {
            double stdx=GLQ.QuadPts[k]; //积分点
            double jx,jy;
            jx=(x[0]+x[1])/2.+(x[1]-x[0])*stdx/2.; jy=(y[0]+y[1])/2.+(y[1]-y[0])*stdx/2.;
            for(i=0;i<polydim2;i++)
            {
                GetBasisFuntion2(E2ID, i, jx, jy, v1);
                for(j=0;j<polydim2;j++)
                {
                    GetBasisFuntion2(E1ID, j, jx, jy, v2);
                    C[i][j]-=(v1[0]*v2[0]+v1[1]*v2[1])*GLQ.Weights[k]/2;
                }
            }
        }
    }
}
*/

void MFEP1NCRS::GetMixedCH1(int ElemID,double **C)
{
    int i,j,k,m,n;
    for(i=0;i<polydim2;i++)
        for(j=0;j<polydim2;j++)
            C[i][j]=0;
    int nv=ms.ElementVertex_Num[ElemID];
    double x[3],y[3],v1[6],v2[6];
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    double xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    for(j=0;j<nv;j++)
    {
        x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
        y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];//积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            for(m=0;m<polydim2;m++)
            {
                GetBasisFuntion2(ElemID, m, jx, jy, v1);
                for(n=0;n<polydim2;n++)
                {
                    GetBasisFuntion2(ElemID, n, jx, jy, v2);
                    C[m][n]+=(v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4]+v1[5]*v2[5])*vol_K*jw;
                }
            }
        }
    }
    delete [] xpos; delete [] ypos;
}

void MFEP1NCRS::GetMixedCL2(int ElemID,double **C)
{
    int i,j,k,m,n;
    for(i=0;i<polydim2;i++)
        for(j=0;j<polydim2;j++)
            C[i][j]=0;
    int nv=ms.ElementVertex_Num[ElemID];
    double x[3],y[3],v1[6],v2[6];
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    double xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    for(j=0;j<nv;j++)
    {
        x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
        y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];//积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            for(m=0;m<polydim2;m++)
            {
                GetBasisFuntion2(ElemID, m, jx, jy, v1);
                for(n=0;n<polydim2;n++)
                {
                    GetBasisFuntion2(ElemID, n, jx, jy, v2);
                    C[m][n]+=(v1[0]*v2[0]+v1[1]*v2[1])*vol_K*jw;
                }
            }
        }
    }
    delete [] xpos; delete [] ypos;
}

void MFEP1NCRS::GetRHS1(int ElemID,FunctionP Source,double *LocF)
{
    int i,j,k,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim1;i++)    LocF[i]=0;
    double fval[2],value[5],x[3],y[3];
    double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<nv;i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];   ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    for(i=0;i<nv;i++)
    {
        x[0]=xE;x[1]=xpos[i];x[2]=xpos[(i+1)%nv];
        y[0]=yE;y[1]=ypos[i];y[2]=ypos[(i+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];     //积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            Source(jx,jy,fval);
            for(j=0;j<polydim1;j++)
            {
                GetBasisFuntion1(ElemID, j, jx, jy, value);
                LocF[j]-=hE*hE*(fval[0]*value[3]+fval[1]*value[4])*vol_K*jw;
            }
        }
    }
    delete []xpos; delete []ypos;
}

void MFEP1NCRS::GetRHS2(int ElemID,FunctionP Source,double *LocF)
{
    int i,j,k,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim2;i++)    LocF[i]=0;
    double fval[2],value[6],x[3],y[3];
    double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置
    double xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<nv;i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];   ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    for(i=0;i<nv;i++)
    {
        x[0]=xE;x[1]=xpos[i];x[2]=xpos[(i+1)%nv];
        y[0]=yE;y[1]=ypos[i];y[2]=ypos[(i+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];     //积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            Source(jx,jy,fval);
            for(j=0;j<polydim2;j++)
            {
                GetBasisFuntion2(ElemID, j, jx, jy, value);
                LocF[j]+=(fval[0]*value[0]+fval[1]*value[1])*vol_K*jw;
            }
        }
    }
    delete []xpos; delete []ypos;
}

void MFEP1NCRS::GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval)
{
    int i,j,k,n=0,en;
    double value[2],ex[2],ey[2];
    
    if(dof2.Num_PerEdge>0)
    {
        GaussLegendreQuadrature GLQ(3);
        double jx,jy,stdx;
        for(i=0;i<ms.B_Edge_Num;i++)
        {
            en=ms.Edge[ms.B_Edge[i] ][0];
            ex[0]=ms.Node[en][0]; ey[0]=ms.Node[en][1];
            en=ms.Edge[ms.B_Edge[i] ][1];
            ex[1]=ms.Node[en][0]; ey[1]=ms.Node[en][1];
            for(j=0;j<dof2.Num_PerEdge;j++)
            {
                Bdof[n]=dof1.Dof_Num+dof2.EdgeD[ ms.B_Edge[i] ] [j];
                for(k=0;k<GLQ.QuadPtsNum;k++)
                {
                    stdx=GLQ.QuadPts[k]; //积分点
                    jx=(ex[0]+ex[1])/2.+(ex[1]-ex[0])*stdx/2.; jy=(ey[0]+ey[1])/2.+(ey[1]-ey[0])*stdx/2.;
                    BFunc(jx,jy,value);
                    Bdofval[n]+=value[j]/2*GLQ.Weights[k];
                }
                n++;
            }
        }
    }
}

void MFEP1NCRS::GetDof1Val(FunctionP sigma,VectorXd &sigmaI)
{
    int i,j,k,nv;
    double value[3],x[3],y[3];
    if(dof1.Num_PerElement>0)
        for(int ElemID=0;ElemID<ms.Element_Num;ElemID++)
        {
            nv=ms.ElementVertex_Num[ElemID];
            double hE=ms.ElementDiameter[ElemID], xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
            double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置
            for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
            {
                xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];    ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
            }

            for(i=0;i<dof1.Num_PerElement;i++)
            {
                sigmaI[ dof1.ElementD[ElemID][i] ]=0;
                for(j=0;j<nv;j++)
                {
                    x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
                    y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
                    double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
                    for(k=0;k<TQ.QuadPtsNum;k++)
                    {
                        double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];//积分点
                        double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
                        double jw=TQ.Weights[k];//积分权重
                        sigma(jx,jy,value);
                        if(i==0)
                            sigmaI[ dof1.ElementD[ElemID][i] ]+=value[0]*vol_K*jw/ms.ElementMeasure[ElemID];
                        else if(i==1)
                            sigmaI[ dof1.ElementD[ElemID][i] ]+=value[0]*(jx-xE)*2*sqrt(2)/hE*vol_K*jw/ms.ElementMeasure[ElemID];
                        else if(i==2)
                            sigmaI[ dof1.ElementD[ElemID][i] ]+=value[1]*vol_K*jw/ms.ElementMeasure[ElemID];
                        else if(i==3)
                            sigmaI[ dof1.ElementD[ElemID][i] ]+=value[2]*vol_K*jw/ms.ElementMeasure[ElemID];
                        else if(i==4)
                            sigmaI[ dof1.ElementD[ElemID][i] ]+=value[2]*(jy-yE)*2*sqrt(2)/hE*vol_K*jw/ms.ElementMeasure[ElemID];
                    }
                }
            }
            delete []xpos; delete[] ypos;
        }
}

void MFEP1NCRS::GetDof2Val(FunctionP u,VectorXd &uI)
{
    int i,j,k,en;
    double ex[2],ey[2];
    double value[6];
    if(dof2.Num_PerEdge>0)
    {
        GaussLegendreQuadrature GLQ(4);
        double jx,jy,stdx;
        for(i=0;i<ms.Edge_Num;i++)
        {
            en=ms.Edge[i][0];
            ex[0]=ms.Node[en][0]; ey[0]=ms.Node[en][1];
            en=ms.Edge[i][1];
            ex[1]=ms.Node[en][0]; ey[1]=ms.Node[en][1];
            
            for(j=0;j<dof2.Num_PerEdge;j++)
            {
                for(k=0;k<GLQ.QuadPtsNum;k++)
                {
                    stdx=GLQ.QuadPts[k]; //积分点
                    jx=(ex[0]+ex[1])/2.+(ex[1]-ex[0])*stdx/2.; jy=(ey[0]+ey[1])/2.+(ey[1]-ey[0])*stdx/2.;
                    u(jx,jy,value);
                    if(j==0)
                        uI[dof2.EdgeD[i][j]]+=value[j]/2*GLQ.Weights[k];
                    else
                        uI[dof2.EdgeD[i][j]]+=value[3]/2*GLQ.Weights[k];
                }
            }
        }
    }
}

void MFEP1NCRS::CheckBasisFunc(int ElemID,double **A)
{
    int i,j,k,m;
    for(i=0;i<8;i++)
        for(j=0;j<8;j++)
            A[i][j]=0;
    int nv=ms.ElementVertex_Num[ElemID];
    double x[3],y[3],v1[6];
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    GaussLobattoQuadrature GLQ(3);
    for(i=0;i<nv;i++)
    {
        x[0]=xpos[i];x[1]=xpos[(i+1)%nv];
        y[0]=ypos[i];y[1]=ypos[(i+1)%nv];
        
        for(k=0;k<GLQ.QuadPtsNum;k++)
        {
            double stdx=GLQ.QuadPts[k]; //积分点
            double jx,jy;
            jx=(x[0]+x[1])/2.+(x[1]-x[0])*stdx/2.; jy=(y[0]+y[1])/2.+(y[1]-y[0])*stdx/2.;
            for(m=0;m<polydim2;m++)
            {
                GetBasisFuntion2(ElemID, m, jx, jy, v1);
                if(m%2==0)
                    A[m][2*i]+=v1[0]*GLQ.Weights[k]/2;
                else
                    A[m][2*i+1]+=v1[1]*GLQ.Weights[k]/2;
            }
        }
    }
    delete [] xpos; delete [] ypos;
}


