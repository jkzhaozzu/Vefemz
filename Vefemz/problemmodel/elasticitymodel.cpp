#include <iostream>
#include <Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers> 

#include "problemmodel.h"
#include "mathfunction.h"
#include "sparsesolve.h"

const double ElasticityProblem::mu=1.;
const double ElasticityProblem::lambda=1000000.;
void ElasticityProblem::Source(double x,double y,double *value)
{
//	value[0]=-PI*PI*(-(3*mu+lambda)*sin(PI*x)*sin(PI*y)+(mu+lambda)*cos(PI*x)*cos(PI*y));
//	value[1]=-PI*PI*(-(3*mu+lambda)*sin(PI*x)*sin(PI*y)+(mu+lambda)*cos(PI*x)*cos(PI*y));

//	value[0]=4*sin(2*PI*y)*cos(2*PI*x)*PI*PI+3*sin(PI*x)*PI*PI*sin(PI*y)/(1+lambda)-lambda*(-sin(PI*x)*PI*PI*sin(PI*y)/(1+lambda)+cos(PI*x)*PI*PI*cos(PI*y)/(1+lambda))+4*sin(2*PI*y)*PI*PI*(-1+cos(2*PI*x))-cos(PI*x)*PI*PI*cos(PI*y)/(1+lambda);
//	value[1]=-4*cos(2*PI*y)*PI*PI*sin(2*PI*x)-cos(PI*x)*PI*PI*cos(PI*y)/(1+lambda)-4*sin(2*PI*x)*PI*PI*(-1+cos(2*PI*y))+3*sin(PI*x)*PI*PI*sin(PI*y)/(1+lambda)-lambda*(-sin(PI*x)*PI*PI*sin(PI*y)/(1+lambda)+cos(PI*x)*PI*PI*cos(PI*y)/(1+lambda));
    
    value[0]=-8*(x+y)*((3*x*y-2)*(x*x+y*y)+5*(x*y-1)*(x*y-1)-2*x*x*y*y);
    value[1]=-8*(x-y)*((3*x*y+2)*(x*x+y*y)-5*(x*y+1)*(x*y+1)+2*x*x*y*y);
}

void ElasticityProblem::DirichletBoundary(double x,double y,double *value)
{
	value[0]=0; value[1]=0;
}
void ElasticityProblem::Solution(double x,double y,double *value)
{
/*	value[0]=sin(PI*x)*sin(PI*y);
	value[1]=PI*cos(PI*x)*sin(PI*y);
	value[2]=PI*cos(PI*y)*sin(PI*x);
	value[3]=sin(PI*x)*sin(PI*y);
	value[4]=PI*cos(PI*x)*sin(PI*y);
	value[5]=PI*cos(PI*y)*sin(PI*x);

	value[0]=sin(2*PI*y)*(-1+cos(2*PI*x))+1/(1+lambda)*sin(PI*x)*sin(PI*y);//u1
	value[1]=-2*sin(2*PI*y)*sin(2*PI*x)*PI+cos(PI*x)*PI*sin(PI*y)/(1+lambda);//u1x
	value[2]=2*cos(2*PI*y)*PI*(-1+cos(2*PI*x))+sin(PI*x)*cos(PI*y)*PI/(1+lambda);//u1y
	value[3]=-sin(2*PI*x)*(-1+cos(2*PI*y))+1/(1+lambda)*sin(PI*x)*sin(PI*y);//u2
	value[4]=-2*cos(2*PI*x)*PI*(-1+cos(2*PI*y))+cos(PI*x)*PI*sin(PI*y)/(1+lambda);//u2x
	value[5]=2*sin(2*PI*y)*sin(2*PI*x)*PI+sin(PI*x)*cos(PI*y)*PI/(1+lambda);//u2y
*/
    value[0]=-4*y*(1-y*y)*(1-x*x)*(1-x*x)-4./(2+lambda)*x*(1-x*x)*(1-y*y)*(1-y*y);
    value[1]=16*y*(1-y*y)*(1-x*x)*x-4*(1-x*x)*(1-y*y)*(1-y*y)/(2+lambda)+8*x*x*(1-y*y)*(1-y*y)/(2+lambda);
    value[2]=-4*(1-y*y)*(1-x*x)*(1-x*x)+8*y*y*(1-x*x)*(1-x*x)+16*x*(1-x*x)*(1-y*y)*y/(2+lambda);
    value[3]=4*x*(1-x*x)*(1-y*y)*(1-y*y)-4./(2+lambda)*y*(1-y*y)*(1-x*x)*(1-x*x);
    value[4]=4*(1-x*x)*(1-y*y)*(1-y*y)-8*x*x*(1-y*y)*(1-y*y)+16*x*(1-x*x)*(1-y*y)*y/(2+lambda);
    value[5]=-16*y*(1-y*y)*(1-x*x)*x-4*(1-y*y)*(1-x*x)*(1-x*x)/(2+lambda)+8*y*y*(1-x*x)*(1-x*x)/(2+lambda);
}

void ElasticityProblem::Stress(double x,double y,double *value)
{
    value[0]=32*y*(1-y*y)*(1-x*x)*x-8*(1-x*x)*(1-y*y)*(1-y*y)/(2+lambda)+16*x*x*(1-y*y)*(1-y*y)/(2+lambda)+lambda*(-4*(1-x*x)*(1-y*y)*(1-y*y)/(2+lambda)+8*x*x*(1-y*y)*(1-y*y)/(2+lambda)-4*(1-y*y)*(1-x*x)*(1-x*x)/(2+lambda)+8*y*y*(1-x*x)*(1-x*x)/(2+lambda));//sigma11
    value[1]=-4*(1-y*y)*(1-x*x)*(1-x*x)+8*y*y*(1-x*x)*(1-x*x)+32*x*(1-x*x)*(1-y*y)*y/(2+lambda)+4*(1-x*x)*(1-y*y)*(1-y*y)-8*x*x*(1-y*y)*(1-y*y);//sigma12
    value[2]=-32*y*(1-y*y)*(1-x*x)*x-8*(1-y*y)*(1-x*x)*(1-x*x)/(2+lambda)+16*y*y*(1-x*x)*(1-x*x)/(2+lambda)+lambda*(-4*(1-x*x)*(1-y*y)*(1-y*y)/(2+lambda)+8*x*x*(1-y*y)*(1-y*y)/(2+lambda)-4*(1-y*y)*(1-x*x)*(1-x*x)/(2+lambda)+8*y*y*(1-x*x)*(1-x*x)/(2+lambda));//sigma22
}

void ElasticityModel::GetRHSL2B(int ElemID,double *LocFB)
{
	int i,j,k,nv=ms.ElementVertex_Num[ElemID];
	for(i=0;i<polydim;i++)
		LocFB[i]=0;
	double *fv=new double[pde.sourcedim]; double value[2];
	double *x=new double[3];	double *y=new double[3];
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
            for(j=0;j<polydim;j++)
            {
			 	pde.Source(jx,jy,fv);
				SMS.GetValue(j,jx, jy, hE, xE, yE,value);
				LocFB[j]+=(fv[0]*value[0]+fv[1]*value[1])*vol_K*jw; 	
            }
        }
	}
	delete [] fv; delete []x; delete []y; delete []xpos; delete []ypos;
}

void ElasticityModel::GetLocRHS(int ElemID,double *LocF)
{
	VE.GetRHS(ElemID,pde.Source,LocF);
/*	int i,j,EdofNum=dof.Total_Num_PerElement[ElemID];
	double *LocFB=new double[polydim];	double **PistarL2=new double*[polydim];	
	for(i=0;i<polydim;i++)	PistarL2[i]=new double[EdofNum];
	
	GetRHSL2B(ElemID,LocFB);	VE.GetPistar_L2(ElemID,PistarL2);

	for(i=0;i<EdofNum;i++)
	{
		LocF[i]=0;
		for(j=0;j<polydim;j++)
			LocF[i]+=LocFB[j]*PistarL2[j][i];
	}

	for(i=0;i<polydim;i++)	delete PistarL2[i];	
	delete []PistarL2; delete []LocFB;
*/
}

void ElasticityModel::Solve()
{
    GetSystem();
    //deal with boundary condition
	int i,k,BdofNum=ms.B_Node_Num*dof.Num_PerNode+ms.B_Edge_Num*dof.Num_PerEdge,IdofNum=dof.Dof_Num-BdofNum;
    int *Bdof=new int[BdofNum];int *Idof=new int[IdofNum];int *IdofID=new int[dof.Dof_Num];;
	double *BdofVal=new double[BdofNum];
	for(i=0;i<BdofNum;i++)
		BdofVal[i]=0;
	VE.GetBdof_BdofVal(pde.DirichletBoundary,Bdof,BdofVal);

    long N=StiffMatrix.nonZeros();
    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(N);
    
    for(i=0;i<BdofNum;i++)
        uh[Bdof[i]]=BdofVal[i];
    
    RHS-=StiffMatrix*uh.pVector;
    
    VectorXd b(IdofNum),x(IdofNum);
    int col=0,row;
    for(i=0;i<dof.Dof_Num;i++) //get the relation between dof and inner dof
    {
        IdofID[i]=-1;
        if(dof.DType[i]==0)
        {
            b[col]=RHS[i];
            Idof[col]=i;
            IdofID[i]=col;
            col++;
        }
    }
    
    for(col=0;col<IdofNum;col++)
    {
        i=Idof[col];
        for(SparseMatrix<double>::InnerIterator it(StiffMatrix,i);it;++it)
        {
            k=(int) it.row();
            if(dof.DType[k]==0)
            {
                row=IdofID[k];
                tripletList.push_back(Eigen::Triplet<double>(row,col,it.value()));
            }
        }
    }
    SparseMatrix<double> A(IdofNum,IdofNum);
    A.setFromTriplets(tripletList.begin(),tripletList.end());
    
    //direct method
    Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    
    // iterative method
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    //    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    x=solver.solve(b);
    for(i=0;i<IdofNum;i++) uh[Idof[i]]=x[i];
    
    delete []Bdof; delete []BdofVal; delete []Idof; delete []IdofID;

}

void ElasticityModel::GetSystem()
{
	int i,j,k,EdofNum,N=0;
    for(i=0;i<ms.Element_Num;i++)
    {
        j=VE.dof.Total_Num_PerElement[i];  N+=j*j;
    }
    double temp;
    std::vector< Eigen::Triplet<double> > tripletList; //for constructing sparsmatrix
    tripletList.reserve(N);
	
	for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
	{
		EdofNum=dof.Total_Num_PerElement[i];
		double **LocA=new double*[EdofNum];double **LocAdiv=new double*[EdofNum];
		double *LocF=new double[EdofNum];
		for(j=0;j<EdofNum;j++)
		{
			LocA[j]=new double[EdofNum];	LocAdiv[j]=new double[EdofNum];
		}
		VE.GetA_H1(i,LocA);		VE.GetA_div(i,LocAdiv);
		GetLocRHS(i,LocF);

		for(j=0;j<EdofNum;j++)//组装
		{
			RHS[ dof.TotalD[i][j] ]+=LocF[j];//组装右端项	
            for(k=0;k<EdofNum;k++)
            {
                temp=pde.mu*LocA[j][k]+(pde.lambda+pde.mu)*LocAdiv[j][k];
                tripletList.push_back(Eigen::Triplet<double>(dof.TotalD[i][j],dof.TotalD[i][k],temp));
            }
		}

		for(j=0;j<EdofNum;j++)
		{
			delete [] LocA[j];delete [] LocAdiv[j];
		}
		delete []LocA; delete [] LocAdiv; delete []LocF;			
	}
    StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());

}

double ElasticityModel::EnergyError()
{
    VectorXd X=uI.pVector-uh.pVector;
    double value=X.transpose()*StiffMatrix*X;
    return sqrt(value);
}

double ElasticityModel::MaxError()
{
	double value=0;
	int i;
	for(i=0;i<dof.Dof_Num;i++)
		if(value<fabs(uI[i]-uh[i]))
			value=fabs(uI[i]-uh[i]);
	return value;
}
