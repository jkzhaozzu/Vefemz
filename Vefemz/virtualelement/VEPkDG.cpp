
#include <iostream>
#include "virtualelement.h"

void VEPkDG::GetB_L2(int ElemID,double **B)
{
	int i,j,EdofNum=dof.Total_Num_PerElement[ElemID];
	for(i=0;i<polydim;i++)
        for(j=0;j<EdofNum;j++)
			if(i==j)
				B[i][j]=ms.ElementMeasure[ElemID];
			else
				B[i][j]=0;
}

void VEPkDG::GetG_L2(int ElemID,double **G)
{
    int i,j,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim;i++)
        for(j=0;j<polydim;j++)
            G[i][j]=0;
    
    double *xpos=new double[nv];	double *ypos=new double[nv];//存储单元位置  
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
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];  //积分点         
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

void VEPkDG::GetPistar_L2(int ElemID,double ** Pistar)//Pistar=G^-1*B
{
	double **B,**G,*X,*R;
	int i,j,EdofNum=dof.Total_Num_PerElement[ElemID];
	X=new double[polydim]; R=new double[polydim];
	B=new double*[polydim]; G=new double*[polydim];
    for(i=0;i<polydim;i++)
	{
		B[i]=new double[EdofNum];	G[i]=new double[polydim];
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

void VEPkDG::GetA_L2(int ElemID,double **A)
{
	int i,j,k,EdofNum=dof.Total_Num_PerElement[ElemID];
	double **B=new double*[polydim];
	double **G=new double*[polydim];
	double **GB=new double*[polydim];
	double **Pistar=new double*[polydim];
	for(i=0;i<polydim;i++)
	{
		B[i]=new double[EdofNum];
		G[i]=new double[polydim];
		GB[i]=new double[EdofNum];
		Pistar[i]=new double[EdofNum];
	}

	GetB_L2(ElemID,B); GetG_L2(ElemID,G); GetPistar_L2(ElemID,Pistar); 
	for(i=0;i<polydim;i++)
		for(j=0;j<EdofNum;j++)
		{
			GB[i][j]=0;
			for(k=0;k<polydim;k++)
				GB[i][j]+=G[i][k]*Pistar[k][j];
		}
	for(i=0;i<EdofNum;i++)
		for(j=0;j<EdofNum;j++)
		{
			A[i][j]=0;
			for(k=0;k<polydim;k++)
				A[i][j]+=Pistar[k][i]*GB[k][j];								
		}	
	
	for(i=0;i<polydim;i++)
	{
		delete[] B[i];	delete[] G[i];	delete[] GB[i];	delete[] Pistar[i];
	}
	delete[] B; delete[] G; delete[] GB; delete[] Pistar;


}

double VEPkDG::GetIntegralValue(int ElemID,double *X)
{
	int i,EdofNum=dof.Total_Num_PerElement[ElemID],nv=ms.ElementVertex_Num[ElemID];
	double **Pistar=new double*[polydim];
	for(i=0;i<polydim;i++)	Pistar[i]=new double[EdofNum];
	GetPistar_L2(ElemID,Pistar);

	double *xpos=new double[nv];	double *ypos=new double[nv];//存储单元位置  
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    
    double x[3],y[3];
    double vi,value=0;
    
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
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];  //积分点         
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;           
            double jw=TQ.Weights[k];//积分权重          
            for(int m=0;m<EdofNum;m++)
				for(int n=0;n<polydim;n++)
				{
					vi=SMS.GetValue(n, jx, jy, hE, xE, yE);
					value+=X[m]*Pistar[n][m]*vi*vol_K*jw;
				}
        }        
    }
    
	for(i=0;i<polydim;i++)
	{
		delete[] Pistar[i];
	}
	delete[] Pistar; delete[] xpos;  delete[] ypos; 
	return value;
}

void VEPkDG::GetDofVal(FunctionP u,double *uI)
{
	int i,j,k,nv;
	double value[3];
	double x[3],y[3];
	
	if(dof.Num_PerElement>0)
		for(int ElemID=0;ElemID<ms.Element_Num;ElemID++)
		{
			nv=ms.ElementVertex_Num[ElemID];			
			double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置			
			double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
			for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
			{
				xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0];	ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
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
						double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];  //积分点             
						double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;						
						double jw=TQ.Weights[k];//积分权重
						u(jx,jy,value);
						uI[dof.ElementD[ElemID][i] ]+=value[0]*SMS.GetValue(i, jx, jy, hE, xE, yE)*vol_K*jw/ms.ElementMeasure[ElemID];
					}
				}
			}
			delete []xpos; delete[] ypos;
		}
}
