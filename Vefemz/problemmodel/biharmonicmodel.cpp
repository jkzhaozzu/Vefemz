#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include "problemmodel.h"
#include "mathfunction.h"
#include "sparsesolve.h"



void BiharmonicProblem::Source(double x,double y,double *value)
{
    *value=4*PI*PI*PI*PI*( (cos(PI*x)*cos(PI*x)-sin(PI*x)*sin(PI*x))*(cos(PI*y)*cos(PI*y)-3*sin(PI*y)*sin(PI*y))+(cos(PI*x)*cos(PI*x)-3*sin(PI*x)*sin(PI*x))*(cos(PI*y)*cos(PI*y)-sin(PI*y)*sin(PI*y)) );
}

void BiharmonicProblem::DirichletBoundary(double x,double y,double *value)
{
    value[0]=0; value[1]=0;
}

void BiharmonicProblem::DirichletBoundary2(double x,double y,double *value)
{
    value[0]=0; value[1]=0; value[2]=0;
}

void BiharmonicProblem::Solution(double x,double y,double *value)
{
    value[0]=sin(PI*x)*sin(PI*x)*sin(PI*y)*sin(PI*y);
    value[1]=2*PI*sin(PI*x)*cos(PI*x)*sin(PI*y)*sin(PI*y);
    value[2]=2*PI*sin(PI*x)*sin(PI*x)*sin(PI*y)*cos(PI*y);
}
void BiharmonicModel::GetRHSL2B(int ElemID,double *LocFB)
{
    int i,j,k,nv=ms.ElementVertex_Num[ElemID],polydim1=VE.p*(VE.p-1)/2;
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
                pde.Source(jx,jy,fv);
                LocFB[j]+=fv[0]*SMS.GetValue(j,jx, jy, hE, xE, yE)*vol_K*jw;
            }
        }
    }
    delete [] fv; delete []x; delete []y; delete []xpos; delete []ypos;
}

void BiharmonicModel::GetLocRHS(int ElemID,double *LocF)
{
    VE.GetRHS(ElemID,pde.Source,LocF);
    /*    int i,j,EdofNum=dof.Total_Num_PerElement[ElemID],polydim1=VE.p*(VE.p-1)/2,nv=ms.ElementVertex_Num[ElemID];
     double *LocFB=new double[polydim1]; double *XX=new double[polydim1];
     double **PistarH1=new double*[polydim];    double **G=new double*[polydim];
     for(i=0;i<polydim;i++)
     {
     PistarH1[i]=new double[EdofNum];G[i]=new double[polydim];
     }
     GetRHSL2B(ElemID,LocFB);    VE.GetPistar_H1(ElemID,PistarH1);VE.GetG_L2(ElemID,G);
     GaussSolve(polydim1,G,LocFB,XX);
     
     double x[3],y[3],temp;
     double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
     double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
     for(i=0;i<ms.ElementVertex_Num[ElemID];i++)
     {
     xpos[i]=ms.Node[ ms.Element[ElemID][i] ][0]; ypos[i]=ms.Node[ ms.Element[ElemID][i] ][1];
     }
     
     for(i=0;i<EdofNum;i++)    LocF[i]=0;
     for(i=0;i<polydim1;i++)
     {
     if(i<dof.Num_PerElement)    LocF[nv+nv*dof.Num_PerEdge+i]+=XX[i]*ms.ElementMeasure[ElemID];
     else
     {
     for(j=0;j<nv;j++)
     {
     x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
     y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
     double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
     for(int k=0;k<TQ.QuadPtsNum;k++)
     {
     double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2]; //积分点
     double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
     double jw=TQ.Weights[k];//积分权重
     temp=XX[i]*SMS.GetValue(i, jx, jy, hE, xE, yE);
     for(int n=0;n<EdofNum;n++)
     for(int m=0;m<polydim;m++)
     LocF[n]+=temp*PistarH1[m][n]*SMS.GetValue(m, jx, jy, hE, xE, yE)*vol_K*jw;
     }
     }
     }
     }
     
     for(i=0;i<polydim;i++)
     {
     delete[] PistarH1[i];    delete[] G[i];
     }
     delete []PistarH1; delete []LocFB; delete[] XX;
     */
}

void BiharmonicModel::Solve()
{
    GetSystem();
    //deal with boundary condition
    int i,k,BdofNum=ms.B_Node_Num*dof.Num_PerNode+ms.B_Edge_Num*dof.Num_PerEdge,IdofNum=dof.Dof_Num-BdofNum;
    int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dof.Dof_Num];
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

void BiharmonicModel::Solve2()
{
    GetSystem();
    //deal with boundary condition
    int i,k,BdofNum=ms.B_Node_Num*dof.Num_PerNode+ms.B_Edge_Num*dof.Num_PerEdge,IdofNum=dof.Dof_Num-BdofNum;
    int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dof.Dof_Num];
    double *BdofVal=new double[BdofNum];
    for(i=0;i<BdofNum;i++)
        BdofVal[i]=0;
    VE.GetBdof_BdofVal(pde.DirichletBoundary2,Bdof,BdofVal);
    
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

void BiharmonicModel::GetSystem()
{
    int i,j,k,EdofNum,N=0;
    for(i=0;i<ms.Element_Num;i++)
    {
        j=VE.dof.Total_Num_PerElement[i];  N+=j*j;
    }
    std::vector< Eigen::Triplet<double> > tripletList; //for constructing sparsmatrix
    tripletList.reserve(N);
    
    for(i=0;i<ms.Element_Num;i++)
    {
        EdofNum=dof.Total_Num_PerElement[i];
        double **LocA=new double*[EdofNum];
        double *LocF=new double[EdofNum];
        for(j=0;j<EdofNum;j++)        LocA[j]=new double[EdofNum];
        VE.GetA_H2(i,LocA);
        GetLocRHS(i,LocF);
        
        for(j=0;j<EdofNum;j++)
        {
            RHS[ dof.TotalD[i][j] ]+=LocF[j];
            for(k=0;k<EdofNum;k++)
                tripletList.push_back(Eigen::Triplet<double>(dof.TotalD[i][j],dof.TotalD[i][k],LocA[j][k]));
            
        }
        for(j=0;j<EdofNum;j++)
            delete [] LocA[j];
        delete []LocA; delete []LocF;
    }
    StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
}

double BiharmonicModel::EnergyError()
{
    VectorXd X=uI.pVector-uh.pVector;
    double value=X.transpose()*StiffMatrix*X;
    return sqrt(value);
}
double BiharmonicModel::MaxError()
{
    double value=0;
    int i;
    for(i=0;i<dof.Dof_Num;i++)
        if(value<fabs(uI[i]-uh[i]))
            value=fabs(uI[i]-uh[i]);
    return value;
}

