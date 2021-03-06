#include <stdio.h>
#include <iostream>
#include "virtualelement.h"
#include "mathfunction.h"

void VEPkNCNodeType::GetB_H1(int ElemID,double **B)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID],EdofNum=dof.Total_Num_PerElement[ElemID];
    
    double e1x[2],e1y[2],n1[2];
    double d1,vx,vy,vxx,vyy;
    
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
        if(i==0)
        {
            if(p==1)
                for(int j=0;j<nv;j++)
				{
                    e1x[0]=xpos[j];	e1x[1]=xpos[(j+1)%nv];
					e1y[0]=ypos[j];	e1y[1]=ypos[(j+1)%nv];
					d1=sqrt((e1x[1]-e1x[0])*(e1x[1]-e1x[0])+(e1y[1]-e1y[0])*(e1y[1]-e1y[0]));			
					B[i][j]+=d1; 
				}
            else
                B[i][nv*dof.Num_PerEdge]=1.;
        }
        else
            for(j=0;j<EdofNum;j++)
            {     
				if(dof.Num_PerEdge>0)
                    if(j>=0&&j<nv*dof.Num_PerEdge)
                    {
                        int en=j/dof.Num_PerEdge;//present edge
                        int ek=j%dof.Num_PerEdge;//present quadrature points
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
                        vx=SMS.GetDerivative(i,1,0,jx,jy,hE,xE,yE);
						vy=SMS.GetDerivative(i,0,1,jx,jy,hE,xE,yE);
                        B[i][j]+=(vx*n1[0]+vy*n1[1])*d1/2*GLQ.Weights[ek];
                    }
				
				if(dof.Num_PerElement>0)
					if(j>=nv*dof.Num_PerEdge) //int_K phi*(-Delta p)dx
					{
						int jp=j-nv*dof.Num_PerEdge;
						int j_row=0,j_col=0,k;
						for(k=0;k<p-1;k++)
							if(jp<(k+2)*(k+1)/2)
							{
								j_row=k;	j_col=jp-(k+1)*k/2;
								break;
							}
						vxx=SMS.GetDerivative(i,j_row-j_col+2,j_col,xE,yE,hE,xE,yE)/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
						vyy=SMS.GetDerivative(i,j_row-j_col,j_col+2,xE,yE,hE,xE,yE)/factorial(j_row-j_col)/factorial(j_col)*ms.ElementMeasure[ElemID];
						for(k=0;k<j_row;k++)
						{
							vxx*=hE;  vyy*=hE;
						}
						B[i][j]-=vxx;	B[i][j]-=vyy;
					}
			}
    }
	delete [] xpos; delete [] ypos;
}

void VEPkNCNodeType::GetD(int ElemID,double **D)
{
    int i,j;  
    for(i=0;i<dof.Total_Num_PerElement[ElemID];i++)//initialize
        for(j=0;j<polydim;j++)
            D[i][j]=0;
    
    int nv=ms.ElementVertex_Num[ElemID];
    double x[3],y[3];
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
                    D[i*dof.Num_PerEdge+k][j]=SMS.GetValue(j,jx,jy,hE,xE,yE);
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
						D[nv*dof.Num_PerEdge+i][m]+=SMS.GetValue(m, jx, jy, hE, xE, yE)*SMS.GetValue(i, jx, jy, hE, xE, yE)*vol_K*jw/ms.ElementMeasure[ElemID];
                    
                    }
                }
            }
    delete [] xpos; delete [] ypos;
}

void VEPkNCNodeType::GetG_H1(int ElemID,double **G)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim;i++)
        for(j=0;j<polydim;j++)
            G[i][j]=0;
    
    double *xpos=new double[nv];//存储单元位置
    double *ypos=new double[nv];
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    
    double x[3],y[3];
    double vix,viy,vjx,vjy;
    
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
				G[0][j]+=SMS.GetValue(j,(x[1]+x[2])/2.,(y[1]+y[2])/2.,hE,xE,yE)*d;
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
                    G[0][j]+=SMS.GetValue(j, jx, jy, hE, xE, yE)*vol_K*jw/ms.ElementMeasure[ElemID];
        
            for(int m=1;m<polydim;m++)
            {
				vix=SMS.GetDerivative(m,1,0, jx, jy, hE, xE, yE);
				viy=SMS.GetDerivative(m,0,1, jx, jy, hE, xE, yE);
                for(int n=1;n<polydim;n++)
                {
                    vjx=SMS.GetDerivative(n,1,0, jx, jy, hE, xE, yE);
					vjy=SMS.GetDerivative(n,0,1, jx, jy, hE, xE, yE);
                    G[m][n]+=(vix*vjx+viy*vjy)*vol_K*jw;
                }				
            }
        }
        
    }
    delete[] xpos;  delete[] ypos; 
}

void VEPkNCNodeType::GetPistar_H1(int ElemID,double ** Pistar)//Pistar=G^-1*B
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

void VEPkNCNodeType::GetPi_H1(int ElemID,double **Pi) //Pi=D*Pistar
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

void VEPkNCNodeType::GetB_L2(int ElemID,double **B)
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
			B[(i+1)*i/2+j][nv*dof.Num_PerEdge+(i+1)*i/2+j]=ms.ElementMeasure[ElemID];
	for(i=p-1;i<p+1;i++)
		for(j=0;j<=i;j++)
			for(m=0;m<EdofNum;m++)
				for(n=0;n<polydim;n++)
					B[(i+1)*i/2+j][m]+=GL2[(i+1)*i/2+j][n]*PistarH1[n][m];
				
	for(i=0;i<polydim;i++)	
	{
		delete[] PistarH1[i];delete[] GL2[i];
	}
	delete[] PistarH1; delete[] GL2; 
}

void VEPkNCNodeType::GetG_L2(int ElemID,double **G)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim;i++)
        for(j=0;j<polydim;j++)
            G[i][j]=0;
    
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

void VEPkNCNodeType::GetPistar_L2(int ElemID,double ** Pistar)//Pistar=G^-1*B
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

void VEPkNCNodeType::GetPi_L2(int ElemID,double **Pi) //Pi=D*Pistar
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

void VEPkNCNodeType::GetA_H1(int ElemID,double **A)
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
			if(i>0)
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

void VEPkNCNodeType::GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval)
{
	int i,j,n=0,en;
	double value[1],ex[2],ey[2];
	
	if(dof.Num_PerEdge>0)
	{
		GaussLegendreQuadrature GLQ(p);
		double jx,jy,stdx;
		for(i=0;i<ms.B_Edge_Num;i++)
		{
			en=ms.Edge[ms.B_Edge[i] ][0];
			ex[0]=ms.Node[en][0]; ey[0]=ms.Node[en][1];
			en=ms.Edge[ms.B_Edge[i] ][1];
			ex[1]=ms.Node[en][0]; ey[1]=ms.Node[en][1];
		
			for(j=0;j<dof.Num_PerEdge;j++)
			{
				Bdof[n]=dof.EdgeD[ ms.B_Edge[i] ] [j];
				stdx=GLQ.QuadPts[j]; //积分点
                jx=(ex[0]+ex[1])/2.+(ex[1]-ex[0])*stdx/2.; jy=(ey[0]+ey[1])/2.+(ey[1]-ey[0])*stdx/2.;
				BFunc(jx,jy,value);
				Bdofval[n]=value[0];
				n++;
			}
		}
	}
}

void VEPkNCNodeType::GetDofVal(FunctionP u,double *uI)
{
	int i,j,k,en,nv;
	double ex[2],ey[2];
	double value[3];
	double x[3],y[3];
	if(dof.Num_PerEdge>0)
	{
		GaussLegendreQuadrature GLQ(p);
		double jx,jy,stdx;
		for(i=0;i<ms.Edge_Num;i++)
		{
			en=ms.Edge[i][0];
			ex[0]=ms.Node[en][0]; ey[0]=ms.Node[en][1];
			en=ms.Edge[i][1];
			ex[1]=ms.Node[en][0]; ey[1]=ms.Node[en][1];
		
			for(j=0;j<dof.Num_PerEdge;j++)
			{
				
				stdx=GLQ.QuadPts[j]; //积分点
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
				xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];
				ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
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
						uI[dof.ElementD[ElemID][i] ]+=value[0]*SMS.GetValue(i, jx, jy, hE, xE, yE)*vol_K*jw/ms.ElementMeasure[ElemID];

					}
				}
			}
			delete []xpos; delete[] ypos;
		}
}

