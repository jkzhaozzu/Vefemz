#include <stdio.h>
#include <iostream>
#include "virtualelement.h"
#include "mathfunction.h"

void VEPkNCVNonPolyDiv::GetB_H1(int ElemID,double **B)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID],EdofNum=dof.Total_Num_PerElement[ElemID];
    
    double e1x[2],e1y[2],n1[2];
    double d1,vx[2],vy[2],vxx[2],vyy[2];
    
    double *xpos=new double[nv];//存储单元位置
    double *ypos=new double[nv];
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];
        ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }

   for(i=0;i<polydim;i++)
        for(j=0;j<EdofNum;j++)
            B[i][j]=0;
    
    GaussLegendreQuadrature GLQ(p);
    int *eflag=ms.ElementEdgeFlag[ElemID];
    for(i=0;i<polydim;i++)
    {
        if(i<2)
        {
            if(p==1)
                for(int j=0;j<nv;j++)
				{
                    e1x[0]=xpos[j];	e1x[1]=xpos[(j+1)%nv];
					e1y[0]=ypos[j];	e1y[1]=ypos[(j+1)%nv];
					d1=sqrt((e1x[1]-e1x[0])*(e1x[1]-e1x[0])+(e1y[1]-e1y[0])*(e1y[1]-e1y[0]));
					if(i==0)
						B[i][2*j]+=d1;
					else
						B[i][2*j+1]+=d1;
				}
            else
			{
				if(i==0)
					B[i][nv*dof.Num_PerEdge]=1.;
				else
					B[i][nv*dof.Num_PerEdge+1]=1.;
			}
        }
        else
            for(j=0;j<EdofNum;j++)
            {     
				if(dof.Num_PerEdge>0)
                    if(j>=0&&j<nv*dof.Num_PerEdge)
                    {
                        int en=j/dof.Num_PerEdge;//present edge
                        int ephi=j%dof.Num_PerEdge;//present basic function on edge en
						int ek=ephi/2;//present quadrature points
                        e1x[0]=xpos[en];    e1x[1]=xpos[(en+1)%nv];
                        e1y[0]=ypos[en];    e1y[1]=ypos[(en+1)%nv];
                        d1=sqrt((e1x[1]-e1x[0])*(e1x[1]-e1x[0])+(e1y[1]-e1y[0])*(e1y[1]-e1y[0]));
                        n1[0]=(e1y[1]-e1y[0])/d1;    n1[1]=(e1x[0]-e1x[1])/d1;
                        
                        double stdx=GLQ.QuadPts[ek]; //积分点
                        double jx,jy;
                        if(eflag[en]==1)
                        {
                            jx=(e1x[0]+e1x[1])/2.+(e1x[1]-e1x[0])*stdx/2.;
                            jy=(e1y[0]+e1y[1])/2.+(e1y[1]-e1y[0])*stdx/2.;
                        }
                        else
                        {
                            jx=(e1x[0]+e1x[1])/2.+(e1x[0]-e1x[1])*stdx/2.;
                            jy=(e1y[0]+e1y[1])/2.+(e1y[0]-e1y[1])*stdx/2.;
                        }
						
                        SMS.GetDerivative(i,1,0,jx,jy,hE,xE,yE,vx);
						SMS.GetDerivative(i,0,1,jx,jy,hE,xE,yE,vy);
						if(ephi%2==0)
							B[i][j]+=(vx[0]*n1[0]+vy[0]*n1[1])*d1/2*GLQ.Weights[ek];
						else
							B[i][j]+=(vx[1]*n1[0]+vy[1]*n1[1])*d1/2*GLQ.Weights[ek];
                    }
				
				if(dof.Num_PerElement>0)
					if(j>=nv*dof.Num_PerEdge) //int_K phi*(-Delta p)dx
					{
						int jp=j-nv*dof.Num_PerEdge;
						int j_row=0,j_col=0,k;
						for(k=0;k<p-1;k++)
							if(jp<(k+2)*(k+1))
							{
								j_row=k;	j_col=(jp-(k+1)*k)/2;
								break;
							}
						SMS.GetDerivative(i,j_row-j_col+2,j_col,xE,yE,hE,xE,yE,vxx);
						vxx[0]=vxx[0]/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
						vxx[1]=vxx[1]/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
						SMS.GetDerivative(i,j_row-j_col,j_col+2,xE,yE,hE,xE,yE,vyy);
						vyy[0]=vyy[0]/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
						vyy[1]=vyy[1]/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
						for(k=0;k<j_row;k++)
						{
							vxx[0]*=hE; vxx[1]*=hE;  vyy[0]*=hE; vyy[1]*=hE;
						}
						if(jp%2==0)
						{
							B[i][j]-=vxx[0];	B[i][j]-=vyy[0];
						}
						else
						{
							B[i][j]-=vxx[1];	B[i][j]-=vyy[1];
						}
					}

			}
    }
	delete [] xpos; delete [] ypos;
}

void VEPkNCVNonPolyDiv::GetD(int ElemID,double **D)
{
    int i,j;  
    for(i=0;i<dof.Total_Num_PerElement[ElemID];i++)//initialize
        for(j=0;j<polydim;j++)
            D[i][j]=0;
    
    int nv=ms.ElementVertex_Num[ElemID];
    double x[3],y[3],v[2],vv[2];
    double *xpos=new double[nv];//存储单元位置
    double *ypos=new double[nv];
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];
        ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
    
    GaussLegendreQuadrature GLQ(p);
    int *eflag=ms.ElementEdgeFlag[ElemID];
    if(dof.Num_PerEdge>0)
        for(i=0;i<nv;i++)
        {
            x[0]=xpos[i];x[1]=xpos[(i+1)%nv];
            y[0]=ypos[i];y[1]=ypos[(i+1)%nv];
            for(int k=0;k<GLQ.QuadPtsNum;k++)
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
                
                for(int j=0;j<polydim;j++)
				{
					SMS.GetValue(j,jx,jy,hE,xE,yE,v);
					D[i*dof.Num_PerEdge+2*k][j]=v[0];
					D[i*dof.Num_PerEdge+2*k+1][j]=v[1];
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
						SMS.GetValue(m, jx, jy, hE, xE, yE,v);
						SMS.GetValue(i, jx, jy, hE, xE, yE,vv);
						D[nv*dof.Num_PerEdge+i][m]+=(v[0]*vv[0]+v[1]*vv[1])*vol_K*jw/ms.ElementMeasure[ElemID];
						
                    
                    }
                }
            }

    delete [] xpos; delete [] ypos;
}

void VEPkNCVNonPolyDiv::GetG_H1(int ElemID,double **G)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim;i++)
        for(j=0;j<polydim;j++)
            G[i][j]=0;
    
    double *xpos=new double[nv];//存储单元位置
    double *ypos=new double[nv];
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    
    double x[3],y[3];
    double vj[2],vix[2],viy[2],vjx[2],vjy[2];
    
    for(i=0;i<nv;i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];
        ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }

    for(i=0;i<nv;i++)
    {
        x[0]=xE;x[1]=xpos[i];x[2]=xpos[(i+1)%nv];
        y[0]=yE;y[1]=ypos[i];y[2]=ypos[(i+1)%nv];
		if(p==1)
		{
			double d=sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
			for(j=0;j<polydim;j++)
			{
				SMS.GetValue(j,(x[1]+x[2])/2.,(y[1]+y[2])/2.,hE,xE,yE,vj);
				G[0][j]+=vj[0]*d;
				G[1][j]+=vj[1]*d;
			}
		}
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(int k=0;k<TQ.QuadPtsNum;k++)
        {
            //积分点
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];           
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            //积分权重
            double jw=TQ.Weights[k];
            if(p>1)
                for(j=0;j<polydim;j++)
				{
					SMS.GetValue(j, jx, jy, hE, xE, yE,vj);
                    G[0][j]+=vj[0]*vol_K*jw/ms.ElementMeasure[ElemID];
					G[1][j]+=vj[1]*vol_K*jw/ms.ElementMeasure[ElemID];
				}
        
            for(int m=2;m<polydim;m++)
            {
				SMS.GetDerivative(m,1,0, jx, jy, hE, xE, yE,vix);
				SMS.GetDerivative(m,0,1, jx, jy, hE, xE, yE,viy);
                for(int n=2;n<polydim;n++)
                {
                    SMS.GetDerivative(n,1,0, jx, jy, hE, xE, yE,vjx);
					SMS.GetDerivative(n,0,1, jx, jy, hE, xE, yE,vjy);
                    G[m][n]+=(vix[0]*vjx[0]+vix[1]*vjx[1]+viy[0]*vjy[0]+viy[1]*vjy[1])*vol_K*jw;
                }				
            }
        }
        
    }
    delete[] xpos;  delete[] ypos; 
}

void VEPkNCVNonPolyDiv::GetPistar_H1(int ElemID,double ** Pistar)//Pistar=G^-1*B
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
	
	GetB_H1(ElemID,B); GetG_H1(ElemID,G);
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
		delete [] B[i];	delete [] G[i];
	}
	delete [] B;	delete [] G;	delete [] X;	delete [] R;	
}

void VEPkNCVNonPolyDiv::GetPi_H1(int ElemID,double **Pi) //Pi=D*Pistar
{
	double **D;
	double **Pistar;
	int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
	D=new double*[EdofNum];	Pistar=new double*[polydim];
	for(i=0;i<EdofNum;i++)
		D[i]=new double[polydim];
	for(i=0;i<polydim;i++)
		Pistar[i]=new double[EdofNum];
	
	GetD(ElemID,D);	GetPistar_H1(ElemID,Pistar);	
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
	delete []D;
	delete []Pistar;
}

void VEPkNCVNonPolyDiv::GetB_L2(int ElemID,double **B)
{
    int i,j,m,n,nv=ms.ElementVertex_Num[ElemID],EdofNum=dof.Total_Num_PerElement[ElemID];
	double **PistarH1=new double*[polydim],**GL2=new double*[polydim];
		
	for(i=0;i<polydim;i++)
	{
		GL2[i]=new double[polydim];
		PistarH1[i]=new double[EdofNum];
	}
	GetPistar_H1(ElemID,PistarH1);
	GetG_L2(ElemID,GL2);

   for(i=0;i<polydim;i++)
        for(j=0;j<EdofNum;j++)
            B[i][j]=0;   
	for(i=0;i<p-1;i++)
		for(j=0;j<=i;j++)
		{
			B[(i+1)*i+2*j][nv*dof.Num_PerEdge+(i+1)*i+2*j]=ms.ElementMeasure[ElemID];
			B[(i+1)*i+2*j+1][nv*dof.Num_PerEdge+(i+1)*i+2*j+1]=ms.ElementMeasure[ElemID];
		}
	for(i=p-1;i<p+1;i++)
		for(j=0;j<=i;j++)
			for(m=0;m<EdofNum;m++)
				for(n=0;n<polydim;n++)
				{
					B[(i+1)*i+2*j][m]+=GL2[(i+1)*i+2*j][n]*PistarH1[n][m];
					B[(i+1)*i+2*j+1][m]+=GL2[(i+1)*i+2*j+1][n]*PistarH1[n][m];
				}
				
	for(i=0;i<polydim;i++)	
	{
		delete[] PistarH1[i];delete[] GL2[i];
	}
	delete[] PistarH1; delete[] GL2; 
}

void VEPkNCVNonPolyDiv::GetG_L2(int ElemID,double **G)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim;i++)
        for(j=0;j<polydim;j++)
            G[i][j]=0;
    
    double *xpos=new double[nv];//存储单元位置
    double *ypos=new double[nv];
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    
    double x[3],y[3];
    double vi[2],vj[2];
    
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];
        ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }

    for(i=0;i<nv;i++)
    {
        x[0]=xE;x[1]=xpos[i];x[2]=xpos[(i+1)%nv];
        y[0]=yE;y[1]=ypos[i];y[2]=ypos[(i+1)%nv];
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
				SMS.GetValue(m, jx, jy, hE, xE, yE,vi);
                for(int n=0;n<polydim;n++)
                {
                    SMS.GetValue(n, jx, jy, hE, xE, yE,vj);
                    G[m][n]+=(vi[0]*vj[0]+vi[1]*vj[1])*vol_K*jw;
                }
            }
        }
        
    }
    delete[] xpos;  delete[] ypos; 
}

void VEPkNCVNonPolyDiv::GetPistar_L2(int ElemID,double ** Pistar)//Pistar=G^-1*B
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
		delete [] B[i];	delete [] G[i];
	}
	delete [] B;	delete [] G;	delete [] X;	delete [] R;	
}

void VEPkNCVNonPolyDiv::GetPi_L2(int ElemID,double **Pi) //Pi=D*Pistar
{
	double **D,**Pistar;
	int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
	D=new double*[EdofNum];	Pistar=new double*[polydim];
	for(i=0;i<EdofNum;i++)
		D[i]=new double[polydim];
	for(i=0;i<polydim;i++)
		Pistar[i]=new double[EdofNum];
	
	GetD(ElemID,D);	GetPistar_L2(ElemID,Pistar);	
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

void VEPkNCVNonPolyDiv::GetA_H1(int ElemID,double **A)
{
	int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
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

	GetB_H1(ElemID,B); GetG_H1(ElemID,G); GetPistar_H1(ElemID,Pistar); GetPi_H1(ElemID,Pi);
	for(i=0;i<polydim;i++)
		for(j=0;j<EdofNum;j++)
		{
			GB[i][j]=0;
			if(i>1)
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
				A[i][j]+=Pi[k][i]*Pi[k][j];					
		}

	
	for(i=0;i<polydim;i++)
	{
		delete[] B[i];	delete[] G[i];	delete[] GB[i];	delete[] Pistar[i];
	}
	for(i=0;i<EdofNum;i++)
		delete [] Pi[i];
	delete[] B; delete[] G; delete[] GB; delete[] Pistar; delete[] Pi;
}

void VEPkNCVNonPolyDiv::GetB_div(int ElemID,double **B)
{
	int i,j,pdim=(p+1)*p/2,nv=ms.ElementVertex_Num[ElemID],EdofNum=dof.Total_Num_PerElement[ElemID];
    
    double ex[2],ey[2],n[2];
    double d,vx,vy;
    
    double *xpos=new double[nv];//存储单元位置
    double *ypos=new double[nv];
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];
        ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }
	for(i=0;i<pdim;i++)
        for(j=0;j<EdofNum;j++)
            B[i][j]=0;
    
    GaussLegendreQuadrature GLQ(p);
	ScaledMonomialSpace sms1d(p-1);
    int *eflag=ms.ElementEdgeFlag[ElemID];

	for(i=0;i<pdim;i++)
	{
		for(j=0;j<EdofNum;j++)
        {     
			if(dof.Num_PerEdge>0)
				if(j>=0&&j<nv*dof.Num_PerEdge)
                {
					int en=j/dof.Num_PerEdge;//present edge
                    int ephi=j%dof.Num_PerEdge;//present basic function on edge en
					int ek=ephi/2;//present quadrature points
                    ex[0]=xpos[en];    ex[1]=xpos[(en+1)%nv];
                    ey[0]=ypos[en];    ey[1]=ypos[(en+1)%nv];
                    d=sqrt((ex[1]-ex[0])*(ex[1]-ex[0])+(ey[1]-ey[0])*(ey[1]-ey[0]));
                    n[0]=(ey[1]-ey[0])/d;    n[1]=(ex[0]-ex[1])/d;
                        
                    double stdx=GLQ.QuadPts[ek]; //积分点
                    double jx,jy;
                    if(eflag[en]==1)
                    {
						jx=(ex[0]+ex[1])/2.+(ex[1]-ex[0])*stdx/2.;
                        jy=(ey[0]+ey[1])/2.+(ey[1]-ey[0])*stdx/2.;
                    }
					else
					{
						jx=(ex[0]+ex[1])/2.+(ex[0]-ex[1])*stdx/2.;
                        jy=(ey[0]+ey[1])/2.+(ey[0]-ey[1])*stdx/2.;
					}
											
					if(ephi%2==0)
						B[i][j]+=n[0]*sms1d.GetValue(i,jx,jy,hE,xE,yE)*d/2*GLQ.Weights[ek];
					else
						B[i][j]+=n[1]*sms1d.GetValue(i,jx,jy,hE,xE,yE)*d/2*GLQ.Weights[ek];
                    }
				
			if(dof.Num_PerElement>0)
				if(j>=nv*dof.Num_PerEdge) //int_K phi*(-Delta p)dx
				{
					int jp=j-nv*dof.Num_PerEdge;
					int j_row=0,j_col=0,k;
					for(k=0;k<p-1;k++)
						if(jp<(k+2)*(k+1))
						{
							j_row=k;	j_col=(jp-(k+1)*k)/2;
							break;
						}
					vx=sms1d.GetDerivative(i,j_row-j_col+1,j_col,xE,yE,hE,xE,yE);
					vx=vx/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
					vy=sms1d.GetDerivative(i,j_row-j_col,j_col+1,xE,yE,hE,xE,yE);
					vy=vy/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
					for(k=0;k<j_row;k++)
					{
						vx*=hE;  vy*=hE; 
					}
					if(jp%2==0)
						B[i][j]-=vx;	
					else
						B[i][j]-=vy;
				
				}
		}
	}
	delete []xpos; delete []ypos;
}

void VEPkNCVNonPolyDiv::GetG_div(int ElemID,double **G)
{
	int i,j,nv=ms.ElementVertex_Num[ElemID],pdim=(p+1)*p/2;
    for(i=0;i<pdim;i++)
        for(j=0;j<pdim;j++)
            G[i][j]=0;
    ScaledMonomialSpace sms1d(p-1);
    double *xpos=new double[nv];//存储单元位置
    double *ypos=new double[nv];
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    
    double x[3],y[3];
    double vi,vj;
    
    for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
    {
        xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];
        ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
    }

    for(i=0;i<nv;i++)
    {
        x[0]=xE;x[1]=xpos[i];x[2]=xpos[(i+1)%nv];
        y[0]=yE;y[1]=ypos[i];y[2]=ypos[(i+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(int k=0;k<TQ.QuadPtsNum;k++)
        {
            //积分点
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];
            
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            //积分权重
            double jw=TQ.Weights[k];
            
            for(int m=0;m<pdim;m++)
            {
				vi=sms1d.GetValue(m, jx, jy, hE, xE, yE);
                for(int n=0;n<pdim;n++)
                {
                    vj=sms1d.GetValue(n, jx, jy, hE, xE, yE);
                    G[m][n]+=vi*vj*vol_K*jw;
                }
            }
        }
        
    }
    delete[] xpos;  delete[] ypos; 
}

void VEPkNCVNonPolyDiv::GetPistar_div(int ElemID,double **Pistar)
{
	double **B,**G,*X,*R;
	int i,j,EdofNum=dof.Total_Num_PerElement[ElemID],pdim=(p+1)*p/2;
	X=new double[pdim]; R=new double[pdim];
	B=new double*[pdim]; G=new double*[pdim];
    for(i=0;i<pdim;i++)
	{
		B[i]=new double[EdofNum];
		G[i]=new double[pdim];
	}
	
	GetB_div(ElemID,B); GetG_div(ElemID,G);
	for(i=0;i<EdofNum;i++)
	{
		for(j=0;j<pdim;j++)
			R[j]=B[j][i];

		GaussSolve(pdim,G,R,X);
		for(j=0;j<pdim;j++)			
			Pistar[j][i]=X[j];
	}
	for(i=0;i<pdim;i++)
	{
		delete [] B[i];	delete [] G[i];
	}
	delete [] B;	delete [] G;	delete [] X;	delete [] R;

}

void VEPkNCVNonPolyDiv::GetA_div(int ElemID,double **A)
{
	int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID],pdim=(p+1)*p/2;
	double **G=new double*[pdim];
	double **GB=new double*[pdim];
	double **Pistar=new double*[pdim];
	for(i=0;i<pdim;i++)
	{
		G[i]=new double[pdim];
		GB[i]=new double[EdofNum];
		Pistar[i]=new double[EdofNum];
	}

	GetG_div(ElemID,G); GetPistar_div(ElemID,Pistar);
	for(i=0;i<pdim;i++)
		for(j=0;j<EdofNum;j++)
		{
			GB[i][j]=0;
			for(k=0;k<pdim;k++)
				GB[i][j]+=G[i][k]*Pistar[k][j];
		}
	for(i=0;i<EdofNum;i++)
		for(j=0;j<EdofNum;j++)
		{
			A[i][j]=0;
			for(k=0;k<pdim;k++)
				A[i][j]+=Pistar[k][i]*GB[k][j];							
		}

	
	for(i=0;i<pdim;i++)
	{
		delete[] G[i];	delete[] GB[i];	delete[] Pistar[i];
	}
	 delete[] G; delete[] GB; delete[] Pistar; 
}


