#include <stdio.h>
#include <iostream>
#include "virtualelement.h"
#include "mathfunction.h"



void VEMorleyType::GetB_H2(int ElemID,double **B)
{
    int i,j,k,nv=ms.ElementVertex_Num[ElemID],EdofNum=dof.Total_Num_PerElement[ElemID];
    double e1x[2],e1y[2],e2x[2],e2y[2],n1[2],n2[2],t1[2],t2[2];
    double d1,d2,vxx,vyy,vxy,vxxx,vxyy,vxxy,vyyy,vx4,vy4,vx2y2;
    
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];  ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    for(i=0;i<polydim;i++)
        for(j=0;j<EdofNum;j++)
            B[i][j]=0;
    
    GaussLobattoQuadrature GLQ(p+1);
    int *eflag=ms.ElementEdgeFlag[ElemID];
    for(i=0;i<polydim;i++)
    {
        if(i==0)
        {
            for(j=0;j<nv;j++)  B[i][j]=1./nv;
        }
        else if(i>0&&i<3)
        {
            for(j=0;j<nv;j++)
            {
                if(j==0)
                {
                    e1x[0]=xpos[nv-1];    e1x[1]=xpos[j];
                    e1y[0]=ypos[nv-1];    e1y[1]=ypos[j];
                }
                else
                {
                    e1x[0]=xpos[j-1];    e1x[1]=xpos[j];
                    e1y[0]=ypos[j-1];    e1y[1]=ypos[j];
                }
                e2x[0]=xpos[j];    e2x[1]=xpos[(j+1)%nv];
                e2y[0]=ypos[j];    e2y[1]=ypos[(j+1)%nv];
                d1=sqrt((e1x[1]-e1x[0])*(e1x[1]-e1x[0])+(e1y[1]-e1y[0])*(e1y[1]-e1y[0]));
                d2=sqrt((e2x[1]-e2x[0])*(e2x[1]-e2x[0])+(e2y[1]-e2y[0])*(e2y[1]-e2y[0]));
                t1[0]=(e1x[1]-e1x[0])/d1;t1[1]=(e1y[1]-e1y[0])/d1;
                t2[0]=(e2x[1]-e2x[0])/d2;t2[1]=(e2y[1]-e2y[0])/d2;
                B[i][j]+=t1[i-1];    B[i][j]-=t2[i-1];
                
                n2[0]=(e2y[1]-e2y[0])/d2;    n2[1]=(e2x[0]-e2x[1])/d2;
                B[i][nv+j*dof.Num_PerEdge+dof.Num_PerEdge/2]=eflag[j]*n2[i-1];
            }
        }
        else
        {
            for(j=0;j<EdofNum;j++)
            {
                if(dof.Num_PerNode>0)
                    if(j<nv)
                    {
                        if(j==0)
                        {
                            e1x[0]=xpos[nv-1];    e1x[1]=xpos[j];
                            e1y[0]=ypos[nv-1];    e1y[1]=ypos[j];
                        }
                        else
                        {
                            e1x[0]=xpos[(j-1)%nv];    e1x[1]=xpos[j];
                            e1y[0]=ypos[(j-1)%nv];    e1y[1]=ypos[j];
                        }
                        e2x[0]=xpos[j];    e2x[1]=xpos[(j+1)%nv];
                        e2y[0]=ypos[j];    e2y[1]=ypos[(j+1)%nv];
                        d1=sqrt((e1x[1]-e1x[0])*(e1x[1]-e1x[0])+(e1y[1]-e1y[0])*(e1y[1]-e1y[0]));
                        d2=sqrt((e2x[1]-e2x[0])*(e2x[1]-e2x[0])+(e2y[1]-e2y[0])*(e2y[1]-e2y[0]));
                        n1[0]=(e1y[1]-e1y[0])/d1;    n1[1]=(e1x[0]-e1x[1])/d1;
                        n2[0]=(e2y[1]-e2y[0])/d2;    n2[1]=(e2x[0]-e2x[1])/d2;
                        t1[0]=(e1x[1]-e1x[0])/d1;t1[1]=(e1y[1]-e1y[0])/d1;
                        t2[0]=(e2x[1]-e2x[0])/d2;t2[1]=(e2y[1]-e2y[0])/d2;
                        
                        vxx=SMS.GetDerivative(i,2,0,e1x[1],e1y[1],hE,xE,yE);
                        vxy=SMS.GetDerivative(i,1,1,e1x[1],e1y[1],hE,xE,yE);
                        vyy=SMS.GetDerivative(i,0,2,e1x[1],e1y[1],hE,xE,yE);
                        B[i][j]+=vxx*t1[0]*n1[0]+vxy*t1[0]*n1[1]+vxy*t1[1]*n1[0]+vyy*t1[1]*n1[1];
                        
                        vxx=SMS.GetDerivative(i,2,0,e2x[0],e2y[0],hE,xE,yE);
                        vxy=SMS.GetDerivative(i,1,1,e2x[0],e2y[0],hE,xE,yE);
                        vyy=SMS.GetDerivative(i,0,2,e2x[0],e2y[0],hE,xE,yE);
                        B[i][j]-=vxx*t2[0]*n2[0]+vxy*t2[0]*n2[1]+vxy*t2[1]*n2[0]+vyy*t2[1]*n2[1];
                    }
                
                if(dof.Num_PerEdge>0)
                    if(j>=nv&&j<nv+nv*dof.Num_PerEdge)
                    {
                        int en=(j-nv)/dof.Num_PerEdge;//present edge
                        int edoftype=(j-nv)%dof.Num_PerEdge;
                        int ek;//present edge dof
                        if(edoftype<dof.Num_PerEdge/2)
                            ek=edoftype;
                        else
                            ek=edoftype-dof.Num_PerEdge/2;
                        e1x[0]=xpos[en];    e1x[1]=xpos[(en+1)%nv];
                        e1y[0]=ypos[en];    e1y[1]=ypos[(en+1)%nv];
                        d1=sqrt((e1x[1]-e1x[0])*(e1x[1]-e1x[0])+(e1y[1]-e1y[0])*(e1y[1]-e1y[0]));
                        n1[0]=(e1y[1]-e1y[0])/d1;    n1[1]=(e1x[0]-e1x[1])/d1;
                        t1[0]=(e1x[1]-e1x[0])/d1;    t1[1]=(e1y[1]-e1y[0])/d1;
                        
                        double jx,jy,tempx,tempy;
                        jx=(e1x[0]+e1x[1])/2.;jy=(e1y[0]+e1y[1])/2.;
                        if(eflag[en]==1)
                        {
                            tempx=(e1x[1]-e1x[0])/2.;tempy=(e1y[1]-e1y[0])/2.;
                        }
                        else
                        {
                            tempx=(e1x[0]-e1x[1])/2.;tempy=(e1y[0]-e1y[1])/2.;
                        }
                        vxx=0;vxy=0;vyy=0;vxxx=0;vxyy=0;vxxy=0;vyyy=0;
                        for(k=0;k<ek+1;k++)
                        {
                            vxx+=SMS.GetDerivative(i,2+k,ek-k,jx,jy,hE,xE,yE)/factorial(k)/factorial(ek-k)*pow(tempx,k)*pow(tempy,ek-k);
                            vxy+=SMS.GetDerivative(i,1+k,1+ek-k,jx,jy,hE,xE,yE)/factorial(k)/factorial(ek-k)*pow(tempx,k)*pow(tempy,ek-k);
                            vyy+=SMS.GetDerivative(i,k,2+ek-k,jx,jy,hE,xE,yE)/factorial(k)/factorial(ek-k)*pow(tempx,k)*pow(tempy,ek-k);
                        }
                        for(k=0;k<ek+1;k++)
                        {
                            vxxx+=SMS.GetDerivative(i,3+k,ek-k,jx,jy,hE,xE,yE)/factorial(k)/factorial(ek-k)*pow(tempx,k)*pow(tempy,ek-k);
                            vxyy+=SMS.GetDerivative(i,1+k,2+ek-k,jx,jy,hE,xE,yE)/factorial(k)/factorial(ek-k)*pow(tempx,k)*pow(tempy,ek-k);
                            vxxy+=SMS.GetDerivative(i,2+k,1+ek-k,jx,jy,hE,xE,yE)/factorial(k)/factorial(ek-k)*pow(tempx,k)*pow(tempy,ek-k);
                            vyyy+=SMS.GetDerivative(i,k,3+ek-k,jx,jy,hE,xE,yE)/factorial(k)/factorial(ek-k)*pow(tempx,k)*pow(tempy,ek-k);
                        }
                        if(edoftype<dof.Num_PerEdge/2)
                        {
                            B[i][j]-=((vxxx+vxyy)*n1[0]+(vxxy+vyyy)*n1[1])*d1;//-int_e v d(\Delta p)/dn ds
                            double temp1=((vxxx*t1[0]+vxxy*t1[1])*n1[0]+(vxxy*t1[0]+vxyy*t1[1])*n1[1])*t1[0];//dp/(dtdndt)
                            double temp2=((vxxy*t1[0]+vxyy*t1[1])*n1[0]+(vxyy*t1[0]+vyyy*t1[1])*n1[1])*t1[1];
                            B[i][j]-=(temp1+temp2)*d1;//-int_e v dp//(dtdndt) ds
                        }
                        else
                            B[i][j]+=eflag[en]*(vxx+vyy-(vxx*t1[0]*t1[0]+2*vxy*t1[0]*t1[1]+vyy*t1[1]*t1[1]));
                    }
                
                if(dof.Num_PerElement>0)
                    if(j>=nv+nv*dof.Num_PerEdge) //int_K v*\Delta^2 pdx
                    {
                        int jp=j-nv-nv*dof.Num_PerEdge;
                        int j_row=0,j_col=0,k=0;
                        for(k=0;k<p-1;k++)
                            if(jp<(k+2)*(k+1)/2)
                            {
                                j_row=k;    j_col=jp-(k+1)*k/2;
                                break;
                            }
                        vx4=SMS.GetDerivative(i,j_row-j_col+4,j_col,xE,yE,hE,xE,yE)/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
                        vx2y2=SMS.GetDerivative(i,j_row-j_col+2,j_col+2,xE,yE,hE,xE,yE)/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
                        vy4=SMS.GetDerivative(i,j_row-j_col,j_col+4,xE,yE,hE,xE,yE)/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
                        for(k=0;k<j_row;k++)
                        {
                            vx4*=hE; vx2y2*=hE; vy4*=hE;
                        }
                        B[i][j]+=vx4+2*vx2y2+vy4;
                    }
            }
        }
    }
    delete [] xpos; delete [] ypos;
}

void VEMorleyType::GetD(int ElemID,double **D)
{
    int i,j,k,m,n;
    for(i=0;i<dof.Total_Num_PerElement[ElemID];i++)//initialize
        for(j=0;j<polydim;j++)
            D[i][j]=0;
    
    int nv=ms.ElementVertex_Num[ElemID];
    double x[3],y[3],temp,vx,vy,d,en[2];
    double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    if(dof.Num_PerNode>0)
        for(i=0;i<nv;i++)
            for(j=0;j<polydim;j++)
                D[i][j]=SMS.GetValue(j,xpos[i],ypos[i],hE,xE,yE);
    
    GaussLobattoQuadrature GLQ(p+1);
    int *eflag=ms.ElementEdgeFlag[ElemID];
    if(dof.Num_PerEdge>0)
        for(i=0;i<nv;i++)
        {
            x[0]=xpos[i];x[1]=xpos[(i+1)%nv];
            y[0]=ypos[i];y[1]=ypos[(i+1)%nv];
            d=sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
            en[0]=(y[1]-y[0])/d;    en[1]=(x[0]-x[1])/d;
            for(k=0;k<GLQ.QuadPtsNum;k++)
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
                
                for(m=0;m<p-1;m++)
                    for(j=0;j<polydim;j++)
                    {
                        temp=1;
                        for(n=0;n<m;n++)  temp*=stdx;
                        if(m<p-2)
                            D[nv+i*dof.Num_PerEdge+m][j]+=temp*SMS.GetValue(j,jx,jy,hE,xE,yE)/2*GLQ.Weights[k];//int_e v t^m ds/|e|
                        vx=SMS.GetDerivative(j,1,0,jx,jy,hE,xE,yE);//int_e grad v\cdot n t^m ds
                        vy=SMS.GetDerivative(j,0,1,jx,jy,hE,xE,yE);
                        D[nv+i*dof.Num_PerEdge+dof.Num_PerEdge/2+m][j]+=eflag[i]*temp*(vx*en[0]+vy*en[1])*d/2*GLQ.Weights[k];
                    }
            }
        }
    
    if(dof.Num_PerElement>0)
        for(i=0;i<dof.Num_PerElement;i++)
            for(j=0;j<nv;j++)
            {
                x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
                y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
                double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
                for(int k=0;k<TQ.QuadPtsNum;k++)
                {
                    //积分点
                    double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];
                    double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
                    //积分权重
                    double jw=TQ.Weights[k];
                    for(int m=0;m<polydim;m++)
                    {
                        D[nv+nv*dof.Num_PerEdge+i][m]+=SMS.GetValue(m, jx, jy, hE, xE, yE)*SMS.GetValue(i, jx, jy, hE, xE, yE)*vol_K*jw/ms.ElementMeasure[ElemID];
                    }
                }
            }
    delete [] xpos; delete [] ypos;
}

void VEMorleyType::GetG_H2(int ElemID,double **G)
{
    int i,j,k,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim;i++)
        for(j=0;j<polydim;j++)
            G[i][j]=0;
    
    double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    
    double x[3],y[3];
    double vixx,vixy,viyy,vjxx,vjxy,vjyy;
    
    for(i=0;i<nv;i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    for(i=0;i<polydim;i++)
    {
        for(j=0;j<ms.ElementVertex_Num[ElemID];j++)
            G[0][i]+=SMS.GetValue(i,xpos[j],ypos[j],hE,xE,yE);
        G[0][i]/=nv;
    }
    GaussLegendreQuadrature GLQ(p);
    for(i=0;i<nv;i++)
    {
        x[0]=xE;x[1]=xpos[i];x[2]=xpos[(i+1)%nv];
        y[0]=yE;y[1]=ypos[i];y[2]=ypos[(i+1)%nv];
        
        double d=sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
        for(k=0;k<GLQ.QuadPtsNum;k++)
        {
            double stdx=GLQ.QuadPts[k]; //积分点
            double jx,jy;
            jx=(x[1]+x[2])/2.+(x[2]-x[1])*stdx/2.; jy=(y[1]+y[2])/2.+(y[2]-y[1])*stdx/2.;
            for(j=0;j<polydim;j++)
            {
                G[1][j]+=SMS.GetDerivative(j,1,0,jx,jy,hE,xE,yE)*d/2*GLQ.Weights[k];
                G[2][j]+=SMS.GetDerivative(j,0,1,jx,jy,hE,xE,yE)*d/2*GLQ.Weights[k];
            }
        }
        
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(int k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];  //积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            for(int m=3;m<polydim;m++)
            {
                vixx=SMS.GetDerivative(m,2,0, jx, jy, hE, xE, yE);
                vixy=SMS.GetDerivative(m,1,1, jx, jy, hE, xE, yE);
                viyy=SMS.GetDerivative(m,0,2, jx, jy, hE, xE, yE);
                for(int n=3;n<polydim;n++)
                {
                    vjxx=SMS.GetDerivative(n,2,0, jx, jy, hE, xE, yE);
                    vjxy=SMS.GetDerivative(n,1,1, jx, jy, hE, xE, yE);
                    vjyy=SMS.GetDerivative(n,0,2, jx, jy, hE, xE, yE);
                    G[m][n]+=(vixx*vjxx+2*vixy*vjxy+viyy*vjyy)*vol_K*jw;
                }
            }
        }
        
    }
    delete[] xpos;  delete[] ypos;
}

void VEMorleyType::GetPistar_H2(int ElemID,double **Pistar)
{
    double **B,**G,*X,*R;
    int i,j,EdofNum=dof.Total_Num_PerElement[ElemID];
    X=new double[polydim]; R=new double[polydim];
    B=new double*[polydim]; G=new double*[polydim];
    for(i=0;i<polydim;i++)
    {
        B[i]=new double[EdofNum];    G[i]=new double[polydim];
    }
    
    GetB_H2(ElemID,B); GetG_H2(ElemID,G);
    for(i=0;i<EdofNum;i++)
    {
        for(j=0;j<polydim;j++)
            R[j]=B[j][i];
        
        GaussSolve(polydim,G,R,X);
        for(j=0;j<polydim;j++)
            Pistar[j][i]=X[j];
    }
    for(i=0;i<polydim;i++)
    {
        delete [] B[i];    delete [] G[i];
    }
    delete [] B;    delete [] G;    delete [] X;    delete [] R;
}

void VEMorleyType::GetPi_H2(int ElemID,double **Pi)
{
    double **D;    double **Pistar;
    int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
    D=new double*[EdofNum];    Pistar=new double*[polydim];
    for(i=0;i<EdofNum;i++)
        D[i]=new double[polydim];
    for(i=0;i<polydim;i++)
        Pistar[i]=new double[EdofNum];
    
    GetD(ElemID,D);    GetPistar_H2(ElemID,Pistar);
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
        {
            Pi[i][j]=0;
            for(k=0;k<polydim;k++)
                Pi[i][j]+=D[i][k]*Pistar[k][j];
        }
    for(i=0;i<EdofNum;i++)
        delete[] D[i];
    for(i=0;i<polydim;i++)
        delete[] Pistar[i];
    delete []D;    delete []Pistar;
}

void VEMorleyType::GetA_H2(int ElemID,double **A)
{
    int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
    double hE=ms.ElementDiameter[ElemID];
    double **B=new double*[polydim];    double **G=new double*[polydim];
    double **GB=new double*[polydim];
    double **Pistar=new double*[polydim];    double **Pi=new double*[EdofNum];
    for(i=0;i<polydim;i++)
    {
        B[i]=new double[EdofNum];        G[i]=new double[polydim];
        GB[i]=new double[EdofNum];        Pistar[i]=new double[EdofNum];
    }
    for(i=0;i<EdofNum;i++)    Pi[i]=new double[EdofNum];
    
    GetB_H2(ElemID,B); GetG_H2(ElemID,G); GetPistar_H2(ElemID,Pistar); GetPi_H2(ElemID,Pi);
    for(i=0;i<polydim;i++)
        for(j=0;j<EdofNum;j++)
        {
            GB[i][j]=0;
            if(i>2)
                for(k=0;k<polydim;k++)
                    GB[i][j]+=G[i][k]*Pistar[k][j];
        }
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
            if(i==j) Pi[i][j]-=1.;
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
        {
            A[i][j]=0;
            for(k=0;k<polydim;k++)
                A[i][j]+=Pistar[k][i]*GB[k][j];
            for(k=0;k<EdofNum;k++)
                A[i][j]+=Pi[k][i]*Pi[k][j]/hE/hE;
        }
    
    for(i=0;i<polydim;i++)
    {
        delete[] B[i];    delete[] G[i];    delete[] GB[i];    delete[] Pistar[i];
    }
    for(i=0;i<EdofNum;i++)
        delete [] Pi[i];
    delete[] B; delete[] G; delete[] GB; delete[] Pistar; delete[] Pi;
}

void VEMorleyType::GetB_L2(int ElemID,double **B)
{
    int i,j,m,n,nv=ms.ElementVertex_Num[ElemID],EdofNum=dof.Total_Num_PerElement[ElemID];
    double **PistarH2=new double*[polydim],**GL2=new double*[polydim];
    
    for(i=0;i<polydim;i++)
    {
        GL2[i]=new double[polydim];
        PistarH2[i]=new double[EdofNum];
    }
    GetPistar_H2(ElemID,PistarH2);
    GetG_L2(ElemID,GL2);
    
    for(i=0;i<polydim;i++)
        for(j=0;j<EdofNum;j++)
            B[i][j]=0;
    for(i=0;i<p-3;i++)
        for(j=0;j<=i;j++)
            B[(i+1)*i/2+j][nv+nv*dof.Num_PerEdge+(i+1)*i/2+j]=ms.ElementMeasure[ElemID];
    for(i=p-3;i<p+1;i++)
        for(j=0;j<=i;j++)
            for(m=0;m<EdofNum;m++)
                for(n=0;n<polydim;n++)
                    B[(i+1)*i/2+j][m]+=GL2[(i+1)*i/2+j][n]*PistarH2[n][m];
    
    for(i=0;i<polydim;i++)
    {
        delete[] PistarH2[i];delete[] GL2[i];
    }
    delete[] PistarH2; delete[] GL2;
}

void VEMorleyType::GetG_L2(int ElemID,double **G)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim;i++)
        for(j=0;j<polydim;j++)
            G[i][j]=0;
    
    
    double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    double x[3],y[3];
    double vi,vj;
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];  ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    for(i=0;i<nv;i++)
    {
        x[0]=xE;x[1]=xpos[i];x[2]=xpos[(i+1)%nv];
        y[0]=yE;y[1]=ypos[i];y[2]=ypos[(i+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(int k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2]; //积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            
            for(int m=0;m<polydim;m++)
            {
                vi=SMS.GetValue(m, jx, jy, hE, xE, yE);
                for(int n=0;n<polydim;n++)
                {
                    vj=SMS.GetValue(n, jx, jy, hE, xE, yE);
                    G[m][n]+=vi*vj*vol_K*jw;
                }
            }
        }
    }
    delete[] xpos;  delete[] ypos;
}

void VEMorleyType::GetPistar_L2(int ElemID,double ** Pistar)//Pistar=G^-1*B
{
    double **B,**G,*X,*R;
    int i,j,EdofNum=dof.Total_Num_PerElement[ElemID];
    X=new double[polydim]; R=new double[polydim];
    B=new double*[polydim]; G=new double*[polydim];
    for(i=0;i<polydim;i++)
    {
        B[i]=new double[EdofNum];
        G[i]=new double[polydim];
    }
    
    GetB_L2(ElemID,B); GetG_L2(ElemID,G);
    for(i=0;i<EdofNum;i++)
    {
        for(j=0;j<polydim;j++)
            R[j]=B[j][i];
        
        GaussSolve(polydim,G,R,X);
        for(j=0;j<polydim;j++)
            Pistar[j][i]=X[j];
    }
    for(i=0;i<polydim;i++)
    {
        delete [] B[i];    delete [] G[i];
    }
    delete [] B;    delete [] G;    delete [] X;    delete [] R;
}

void VEMorleyType::GetPi_L2(int ElemID,double **Pi) //Pi=D*Pistar
{
    double **D,**Pistar;
    int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
    D=new double*[EdofNum];    Pistar=new double*[polydim];
    for(i=0;i<EdofNum;i++)
        D[i]=new double[polydim];
    for(i=0;i<polydim;i++)
        Pistar[i]=new double[EdofNum];
    
    GetD(ElemID,D);    GetPistar_L2(ElemID,Pistar);
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
        {
            Pi[i][j]=0;
            for(k=0;k<polydim;k++)
                Pi[i][j]+=D[i][k]*Pistar[k][j];
        }
    for(i=0;i<EdofNum;i++)
        delete [] D[i];
    for(i=0;i<polydim;i++)
        delete [] Pistar[i];
    delete [] D;
    delete [] Pistar;
}

void VEMorleyType::GetA_L2(int ElemID,double **A)
{
    int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
    double hE=ms.ElementDiameter[ElemID];
    double **B=new double*[polydim];
    double **G=new double*[polydim];
    double **GB=new double*[polydim];
    double **Pistar=new double*[polydim];
    double **Pi=new double*[EdofNum];
    for(i=0;i<polydim;i++)
    {
        B[i]=new double[EdofNum];
        G[i]=new double[polydim];
        GB[i]=new double[EdofNum];
        Pistar[i]=new double[EdofNum];
    }
    for(i=0;i<EdofNum;i++)
        Pi[i]=new double[EdofNum];
    
    GetB_L2(ElemID,B); GetG_L2(ElemID,G); GetPistar_L2(ElemID,Pistar); GetPi_L2(ElemID,Pi);
    for(i=0;i<polydim;i++)
        for(j=0;j<EdofNum;j++)
        {
            GB[i][j]=0;
            for(k=0;k<polydim;k++)
                GB[i][j]+=G[i][k]*Pistar[k][j];
        }
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
            if(i==j) Pi[i][j]-=1.;
    for(i=0;i<EdofNum;i++)
        for(j=0;j<EdofNum;j++)
        {
            A[i][j]=0;
            for(k=0;k<polydim;k++)
                A[i][j]+=Pistar[k][i]*GB[k][j];
            for(k=0;k<EdofNum;k++)
                A[i][j]+=Pi[k][i]*Pi[k][j]*hE*hE;
        }
    
    
    for(i=0;i<polydim;i++)
    {
        delete[] B[i];    delete[] G[i];    delete[] GB[i];    delete[] Pistar[i];
    }
    for(i=0;i<EdofNum;i++)
        delete [] Pi[i];
    delete[] B; delete[] G; delete[] GB; delete[] Pistar; delete[] Pi;
}

void VEMorleyType::GetRHSL2B(int ElemID,FunctionP Source,double *LocFB)
{
    int i,j,k,nv=ms.ElementVertex_Num[ElemID],polydim1;
    if(p==2||p==3)    polydim1=p*(p-1)/2;
    else if(p==4)    polydim1=3;
    else polydim1=(p-2)*(p-3)/2;
    for(i=0;i<polydim1;i++)    LocFB[i]=0;
    double *fv=new double[1];
    double *x=new double[3];    double *y=new double[3];
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
            for(j=0;j<polydim1;j++)
            {
                Source(jx,jy,fv);
                LocFB[j]+=fv[0]*SMS.GetValue(j,jx, jy, hE, xE, yE)*vol_K*jw;
            }
        }
    }
    delete [] fv; delete []x; delete []y; delete []xpos; delete []ypos;
}

void VEMorleyType::GetRHS(int ElemID,FunctionP Source,double *LocF)
{
    int i,j,m,n,EdofNum=dof.Total_Num_PerElement[ElemID],polydim1,nv=ms.ElementVertex_Num[ElemID];
    if(p==2||p==3)    polydim1=p*(p-1)/2;
    else if(p==4)    polydim1=3;
    else polydim1=(p-2)*(p-3)/2;
    
    double *LocFB=new double[polydim1]; double *XX=new double[polydim1];    double **G=new double*[polydim];
    double d1,d2,e1x[2],e1y[2],e2x[2],e2y[2],n2[2],t1[2],t2[2],dsum=0;
    for(i=0;i<polydim;i++)    G[i]=new double[polydim];
    
    GetRHSL2B(ElemID,Source,LocFB);    GetG_L2(ElemID,G);
    GaussSolve(polydim1,G,LocFB,XX);
    
    double x[3],y[3],temp;
    int *eflag=ms.ElementEdgeFlag[ElemID];
    double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<nv;i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    for(j=0;j<nv;j++)
    {
        e2x[0]=xpos[j];    e2x[1]=xpos[(j+1)%nv];
        e2y[0]=ypos[j];    e2y[1]=ypos[(j+1)%nv];
        dsum+=sqrt((e2x[1]-e2x[0])*(e2x[1]-e2x[0])+(e2y[1]-e2y[0])*(e2y[1]-e2y[0]));
    }
    
    for(i=0;i<EdofNum;i++)    LocF[i]=0;
    if(p==2||p==3||p==4)
        for(n=0;n<nv;n++)
        {
            x[0]=xE;x[1]=xpos[n];x[2]=xpos[(n+1)%nv];
            y[0]=yE;y[1]=ypos[n];y[2]=ypos[(n+1)%nv];
            double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
            for(int k=0;k<TQ.QuadPtsNum;k++)
            {
                double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2]; //积分点
                double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
                double jw=TQ.Weights[k];//积分权重
                temp=0;
                for(m=0;m<polydim1;m++)    temp+=XX[m]*SMS.GetValue(m, jx, jy, hE, xE, yE);
                if(p==2||p==3)
                    for(m=0;m<nv;m++)    LocF[m]+=temp/nv*vol_K*jw;
                else
                    LocF[nv+nv*dof.Num_PerEdge]+=temp*vol_K*jw;
                for(j=0;j<nv;j++)
                {
                    if(j==0)
                    {
                        e1x[0]=xpos[nv-1];    e1x[1]=xpos[j];
                        e1y[0]=ypos[nv-1];    e1y[1]=ypos[j];
                    }
                    else
                    {
                        e1x[0]=xpos[j-1];    e1x[1]=xpos[j];
                        e1y[0]=ypos[j-1];    e1y[1]=ypos[j];
                    }
                    e2x[0]=xpos[j];    e2x[1]=xpos[(j+1)%nv];
                    e2y[0]=ypos[j];    e2y[1]=ypos[(j+1)%nv];
                    d1=sqrt((e1x[1]-e1x[0])*(e1x[1]-e1x[0])+(e1y[1]-e1y[0])*(e1y[1]-e1y[0]));
                    d2=sqrt((e2x[1]-e2x[0])*(e2x[1]-e2x[0])+(e2y[1]-e2y[0])*(e2y[1]-e2y[0]));
                    t1[0]=(e1x[1]-e1x[0])/d1;t1[1]=(e1y[1]-e1y[0])/d1;
                    t2[0]=(e2x[1]-e2x[0])/d2;t2[1]=(e2y[1]-e2y[0])/d2;
                    n2[0]=(e2y[1]-e2y[0])/d2;    n2[1]=(e2x[0]-e2x[1])/d2;
                    if(p==2||p==3)
                    {
                        LocF[j]+=temp*(p-2)*((jx-xE)*(t1[0]-t2[0])+(jy-yE)*(t1[1]-t2[1]))/dsum*vol_K*jw;
                        LocF[nv+j*dof.Num_PerEdge+dof.Num_PerEdge/2]=eflag[j]*temp*(p-2)*((jx-xE)*n2[0]+(jy-yE)*n2[1])/dsum*vol_K*jw;
                    }
                    else
                    {
                        LocF[j]+=temp*((jx-xE)*(t1[0]-t2[0])+(jy-yE)*(t1[1]-t2[1]))/dsum*vol_K*jw;
                        LocF[nv+j*dof.Num_PerEdge+dof.Num_PerEdge/2]=eflag[j]*temp*((jx-xE)*n2[0]+(jy-yE)*n2[1])/dsum*vol_K*jw;
                    }
                }
            }
        }
    
    if(p>4)
        for(i=0;i<polydim1;i++)
        {
            LocF[nv+nv*dof.Num_PerEdge+i]+=XX[i]*ms.ElementMeasure[ElemID];
        }
    
    for(i=0;i<polydim;i++)    delete[] G[i];
    delete []LocFB; delete[] XX;
}

void VEMorleyType::GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval)
{
    int i,j,k,m,n=0,en;
    double d,temp,value[2],x,y,ex[2],ey[2];
    
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
        GaussLegendreQuadrature GLQ(p+2);
        double jx,jy,stdx;
        for(i=0;i<ms.B_Edge_Num;i++)
        {
            en=ms.Edge[ms.B_Edge[i] ][0];
            ex[0]=ms.Node[en][0]; ey[0]=ms.Node[en][1];
            en=ms.Edge[ms.B_Edge[i] ][1];
            ex[1]=ms.Node[en][0]; ey[1]=ms.Node[en][1];
            d=sqrt((ex[1]-ex[0])*(ex[1]-ex[0])+(ey[1]-ey[0])*(ey[1]-ey[0]));
            for(j=0;j<dof.Num_PerEdge;j++)
            {
                Bdof[n]=dof.EdgeD[ ms.B_Edge[i] ] [j];
                for(k=0;k<GLQ.QuadPtsNum;k++)
                {
                    stdx=GLQ.QuadPts[k]; //积分点
                    jx=(ex[0]+ex[1])/2.+(ex[1]-ex[0])*stdx/2.; jy=(ey[0]+ey[1])/2.+(ey[1]-ey[0])*stdx/2.;
                    BFunc(jx,jy,value);
                    temp=1;
                    if(j<dof.Num_PerEdge/2)
                    {
                        for(m=0;m<j;m++)    temp*=stdx;
                        Bdofval[n]+=temp*value[0]/2*GLQ.Weights[k];
                    }
                    else
                    {
                        for(m=0;m<j-(dof.Num_PerEdge/2);m++)    temp*=stdx;
                        Bdofval[n]+=temp*value[1]*d/2*GLQ.Weights[k];
                    }
                }
                n++;
            }
        }
    }
}

void VEMorleyType::GetDofVal(FunctionP u,double *uI)
{
    int i,j,k,m,en,nv;
    double px,py,temp,d,ex[2],ey[2],n1[2];
    double value[3];
    double x[3],y[3];
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
        GaussLegendreQuadrature GLQ(p+2);
        double jx,jy,stdx;
        for(i=0;i<ms.Edge_Num;i++)
        {
            en=ms.Edge[i][0];
            ex[0]=ms.Node[en][0]; ey[0]=ms.Node[en][1];
            en=ms.Edge[i][1];
            ex[1]=ms.Node[en][0]; ey[1]=ms.Node[en][1];
            d=sqrt((ex[1]-ex[0])*(ex[1]-ex[0])+(ey[1]-ey[0])*(ey[1]-ey[0]));
            n1[0]=(ey[1]-ey[0])/d;    n1[1]=(ex[0]-ex[1])/d;
            
            for(j=0;j<dof.Num_PerEdge;j++)
            {
                int ek;
                if(j<dof.Num_PerEdge/2)
                    ek=j;
                else
                    ek=j-dof.Num_PerEdge/2;
                for(k=0;k<GLQ.QuadPtsNum;k++)
                {
                    stdx=GLQ.QuadPts[k]; //积分点
                    jx=(ex[0]+ex[1])/2.+(ex[1]-ex[0])*stdx/2.; jy=(ey[0]+ey[1])/2.+(ey[1]-ey[0])*stdx/2.;
                    u(jx,jy,value);
                    temp=1;
                    for(m=0;m<ek;m++)
                        temp*=stdx;
                    if(j<dof.Num_PerEdge/2)
                        uI[dof.EdgeD[i][j]]+=temp*value[0]/2*GLQ.Weights[k];
                    else
                        uI[dof.EdgeD[i][j]]+=temp*(value[1]*n1[0]+value[2]*n1[1])*d/2*GLQ.Weights[k];
                }
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
                for(j=0;j<nv;j++)
                {
                    x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
                    y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
                    double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
                    for(k=0;k<TQ.QuadPtsNum;k++)
                    {
                        //积分点
                        double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];
                        double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
                        //积分权重
                        double jw=TQ.Weights[k];
                        u(jx,jy,value);
                        uI[ dof.ElementD[ElemID][i] ]+=value[0]*SMS.GetValue(i, jx, jy, hE, xE, yE)*vol_K*jw/ms.ElementMeasure[ElemID];
                        
                    }
                }
            }
            delete []xpos; delete[] ypos;
        }
}

