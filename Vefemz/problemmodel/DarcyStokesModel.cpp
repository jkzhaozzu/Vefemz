//
//  DarcyStokesModel.cpp
//  Vefemz
//
//  Created by 张蓓 on 2020/3/12.
//  Copyright © 2020年 Jikun Zhao. All rights reserved.
//

#include <iostream>

#include "problemmodel.h"
#include "mathfunction.h"

const double DarcyStokesProblem::epsilon=0; //pow(2.,-14);
/*
 H1-nonconforming and H(div)-conforming divergence-free element
 光滑解时epsilon=1/4，1/16,1/64 on rectmesh, epsilon=1/4，1/16, 1/64 on poly mesh 稳定项中少个epsilon
 */
void DarcyStokesProblem::Source1(double x,double y,double *value)
{
    value[0]=-cos(2*PI*x)*sin(2*PI*y)+sin(2*PI*y)+epsilon*epsilon*(-8*PI*PI*cos(2*PI*x)*sin(2*PI*y)+4*PI*PI*sin(2*PI*y))+y*y;
    value[1]=sin(2*PI*x)*cos(2*PI*y)-sin(2*PI*x)+epsilon*epsilon*(8*PI*PI*cos(2*PI*y)*sin(2*PI*x)-4*PI*PI*sin(2*PI*x))+2*x*y;

//    value[0]=x*exp(-(x*y)/epsilon)-(-2.*epsilon*y+x*y*y+x*x*x)*exp(-(x*y)/epsilon)-exp(-x/epsilon);
//    value[1]=-y*exp(-(x*y)/epsilon)-(-y*y*y+2*epsilon*x-x*x*y)*exp(-(x*y)/epsilon);
}

void DarcyStokesProblem::Source2(double x,double y,double *value)
{
    value[0]=0;
}

void DarcyStokesProblem::DirichletBoundary(double x,double y,double *value)
{
    value[0]=0; value[1]=0;
//    value[0]=x*exp(-(x*y)/epsilon);
 //   value[1]=-y*exp(-(x*y)/epsilon);
}
void DarcyStokesProblem::Velocity(double x,double y,double *value)
{
    value[0]=-cos(2*PI*x)*sin(2*PI*y)+sin(2*PI*y);
    value[1]=0; value[2]=0;
    value[3]=sin(2*PI*x)*cos(2*PI*y)-sin(2*PI*x);
    value[4]=0; value[5]=0;

 /*   value[0]=x*exp(-(x*y)/epsilon);
    value[1]=0; value[2]=0;
    value[3]=-y*exp(-(x*y)/epsilon);
    value[4]=0; value[5]=0;
*/

}

void DarcyStokesProblem::Pressure(double x,double y,double *value)
{
    value[0]=x*y*y-1./6;    value[1]=0;    value[2]=0;
    
//    value[0]=epsilon*exp(-x/epsilon)-epsilon*epsilon*(1.-exp(-1./epsilon));
//    value[1]=0; value[2]=0;
 
}

void DarcyStokesModel::GetRHSL2B(int ElemID,double *LocFB)
{
    int i,j,k,nv=ms.ElementVertex_Num[ElemID];
    for(i=0;i<polydim;i++)
        LocFB[i]=0;
    double *fv=new double[pde.sourcedim]; double value[2];
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
            for(j=0;j<polydim;j++)
            {
                pde.Source1(jx,jy,fv);
                SMSV.GetValue(j,jx, jy, hE, xE, yE,value);
                LocFB[j]+=(fv[0]*value[0]+fv[1]*value[1])*vol_K*jw;
            }
        }
    }
    delete [] fv; delete []x; delete []y; delete []xpos; delete []ypos;
}

void DarcyStokesModel::GetLocRHS(int ElemID,double *LocF)
{
 //   VE1.GetRHS(ElemID,pde.Source1,LocF);
    int i,j,EdofNum=dof1.Total_Num_PerElement[ElemID];
    double *LocFB=new double[polydim];    double **PistarL2=new double*[polydim];
    for(i=0;i<polydim;i++)    PistarL2[i]=new double[EdofNum];
     
    GetRHSL2B(ElemID,LocFB);    VE1.GetPistar_L2(ElemID,PistarL2);
     
    for(i=0;i<EdofNum;i++)
    {
        LocF[i]=0;
        for(j=0;j<polydim;j++)
            LocF[i]+=LocFB[j]*PistarL2[j][i];
    }
     
    for(i=0;i<polydim;i++)    delete PistarL2[i];
    delete []PistarL2; delete []LocFB;
}

void DarcyStokesModel::Solve()
{
    GetSystem();
    //deal with boundary condition for u
    int i,j,k,dofNum=dof1.Dof_Num+dof2.Dof_Num,BdofNum=ms.B_Node_Num*dof1.Num_PerNode+ms.B_Edge_Num*dof1.Num_PerEdge,
    IdofNum=dofNum-BdofNum;
    int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dofNum];
    double *BdofVal=new double[BdofNum];
    for(i=0;i<BdofNum;i++)
        BdofVal[i]=0;
    VE1.GetBdof_BdofVal(pde.DirichletBoundary,Bdof,BdofVal);
    
    long N=StiffMatrix.nonZeros();
    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(N);
    
    VectorXd x0(dof1.Dof_Num+dof2.Dof_Num);
    x0.setZero();
    for(i=0;i<BdofNum;i++)
    {
        uh[Bdof[i]]=BdofVal[i]; x0(Bdof[i])=BdofVal[i];
    }
    RHS-=StiffMatrix*x0;
    
    VectorXd b(IdofNum),x(IdofNum);
    int col=0,row;
    for(i=0;i<dofNum;i++) //get the relation between dof and inner dof
        if(i<dof1.Dof_Num)
        {
            IdofID[i]=-1;
            if(dof1.DType[i]==0)
            {
                b[col]=RHS[i];  Idof[col]=i; IdofID[i]=col;  col++;
            }
        }
        else
        {
            b[col]=RHS[i];  Idof[col]=i; IdofID[i]=col;  col++;
        }
    
    for(col=0;col<IdofNum;col++)
    {
        i=Idof[col];
        for(SparseMatrix<double>::InnerIterator it(StiffMatrix,i);it;++it)
        {
            k=(int) it.row();
            if(k<dof1.Dof_Num)
            {
                if(dof1.DType[k]==0)
                {
                    row=IdofID[k];
                    tripletList.push_back(Eigen::Triplet<double>(row,col,it.value()));
                }
            }
            else
            {
                row=IdofID[k];
                tripletList.push_back(Eigen::Triplet<double>(row,col,it.value()));
            }
        }
    }
    SparseMatrix<double> A(IdofNum,IdofNum);
    A.setFromTriplets(tripletList.begin(),tripletList.end());
    
    //direct method
    //Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
   // Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    
    // iterative method
   // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    x=solver.solve(b);
    
    for(i=0;i<IdofNum;i++)
    {
        k=Idof[i];
        if(k<dof1.Dof_Num)
            uh[k]=x[i];
        else
            ph[k-dof1.Dof_Num]=x[i];
    }
    
    //add L_0^2 condition for ph
    double mean=ph.GetAverageValue();
    for(i=0;i<ph.ms.Element_Num;i++)
        ph[i*ph.dof.Num_PerElement]-=mean;
    
    delete []Bdof; delete []BdofVal; delete []Idof; delete []IdofID;
}

void DarcyStokesModel::UzawaSolve()
{
    GetSystem1();
    //deal with boundary condition for u
    int i,j,k,dofNum=dof1.Dof_Num,BdofNum=ms.B_Node_Num*dof1.Num_PerNode+ms.B_Edge_Num*dof1.Num_PerEdge,
    IdofNum=dof1.Dof_Num-BdofNum;;
    int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dofNum];
    double *BdofVal=new double[BdofNum];
    for(i=0;i<BdofNum;i++)
        BdofVal[i]=0;
    VE1.GetBdof_BdofVal(pde.DirichletBoundary,Bdof,BdofVal);
 
    long N=StiffMatrixA.nonZeros();
    std::vector< Eigen::Triplet<double> > tripletListA;
    tripletListA.reserve(N);
    N=SMB.nonZeros();
    std::vector< Eigen::Triplet<double> > tripletListB;
    tripletListA.reserve(N);
   
    for(i=0;i<BdofNum;i++)
        uh[Bdof[i]]=BdofVal[i];
    RHS1-=StiffMatrixA*uh.pVector;
    
    VectorXd b(IdofNum);
    int col=0,row;
    for(i=0;i<dofNum;i++) //get the relation between dof and inner dof
    {
        IdofID[i]=-1;
        if(dof1.DType[i]==0)
        {
            b[col]=RHS1[i];  Idof[col]=i; IdofID[i]=col;  col++;
        }
    }
    
    for(col=0;col<IdofNum;col++)
    {
        i=Idof[col];
        for(SparseMatrix<double>::InnerIterator it(StiffMatrixA,i);it;++it)
        {
            k=(int) it.row();
            if(dof1.DType[k]==0)
            {
                row=IdofID[k];
                tripletListA.push_back(Eigen::Triplet<double>(row,col,it.value()));
            }
        }
    }
    for(col=0;col<dof2.Dof_Num;col++)
        for(SparseMatrix<double>::InnerIterator it(SMB,col);it;++it)
        {
            k=(int) it.row();
            if(dof1.DType[k]==0)
            {
                row=IdofID[k];
                tripletListB.push_back(Eigen::Triplet<double>(row,col,it.value()));
            }
        }
    SparseMatrix<double> A(IdofNum,IdofNum),B(IdofNum,dof2.Dof_Num);
    A.setFromTriplets(tripletListA.begin(),tripletListA.end());
    B.setFromTriplets(tripletListB.begin(),tripletListB.end());

 
    //direct method
    //Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    
    // iterative method
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
     Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver2;
    solver2.compute(StiffMatrixC);
    VectorXd x1(IdofNum),x0(IdofNum),y1(dof2.Dof_Num),y0(dof2.Dof_Num);
    x0.setZero();y0.setZero();
    x1.setZero();y1.setZero();
    double err;
    for(i=0;i<100;i++)
    {
        x0=x1;y0=y1;
        x1=solver.solve(b-B*y0);
        y1=y0+solver2.solve(B.transpose()*x1);
        err=(x1-x0).dot(x1-x0)+(y1-y0).dot(y1-y0);
        if(err<0.0000001) break;
    }
 
 //   std::cout<<y0<<std::endl;
    for(i=0;i<IdofNum;i++)
    {
        k=Idof[i]; uh[k]=x1[i];
    }
    for(i=0;i<dof2.Dof_Num;i++)
        ph[i]=y1[i];
    
    //add L_0^2 condition for ph
    double mean=ph.GetAverageValue();
    for(i=0;i<ph.ms.Element_Num;i++)
        ph[i*ph.dof.Num_PerElement]-=mean;
    delete []Bdof; delete []BdofVal; delete []Idof; delete []IdofID;
}

void DarcyStokesModel::GetSystem()
{
    int i,j,k,EdofNum1,EdofNum2,N=0,N1=0,N2=0;
    for(i=0;i<ms.Element_Num;i++)
    {
        j=VE1.dof.Total_Num_PerElement[i]; k=VE2.dof.Total_Num_PerElement[i];
        N+=j*j*2+2*j*k; N1+=j*j*2; N2+=k*k;
    }
    std::vector< Eigen::Triplet<double> > tripletList; //for constructing sparsmatrix
    tripletList.reserve(N);
    std::vector< Eigen::Triplet<double> > tripletList1; //for constructing sparsmatrix A
    tripletList1.reserve(N1);
    std::vector< Eigen::Triplet<double> > tripletList2; //for constructing sparsmatrix C
    tripletList2.reserve(N2);
    for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
    {
        EdofNum1=dof1.Total_Num_PerElement[i];EdofNum2=dof2.Total_Num_PerElement[i];
        double **LocA=new double*[EdofNum1];double **LocAL2=new double*[EdofNum1];
        double **LocAdiv=new double*[EdofNum1];
        double *LocF=new double[EdofNum1];double **LocC=new double*[EdofNum2];
        for(j=0;j<EdofNum1;j++)
        {
            LocA[j]=new double[EdofNum1];   LocAL2[j]=new double[EdofNum1];
            LocAdiv[j]=new double[EdofNum2];
        }
        for(j=0;j<EdofNum2;j++)    LocC[j]=new double[EdofNum2];
        MVE.GetMixedA(i,LocA);        MVE.GetMixedB(i,LocAdiv);
        VE1.GetA_L2(i, LocAL2);
        GetLocRHS(i,LocF);VE2.GetA_L2(i,LocC);
        for(j=0;j<EdofNum1;j++)//组装 A,B
        {
            RHS[ dof1.TotalD[i][j] ]+=LocF[j];//组装右端项
            for(k=0;k<EdofNum1;k++)
            {
                tripletList.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.TotalD[i][k],pde.epsilon*pde.epsilon*LocA[j][k]+LocAL2[j][k]));
                tripletList1.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.TotalD[i][k],pde.epsilon*pde.epsilon*LocA[j][k]+LocAL2[j][k]));
            }
            for(k=0;k<EdofNum2;k++)
            {
                tripletList.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.Dof_Num+dof2.TotalD[i][k],LocAdiv[j][k]));
                tripletList.push_back(Eigen::Triplet<double>(dof1.Dof_Num+dof2.TotalD[i][k],dof1.TotalD[i][j],LocAdiv[j][k]));
            }
        }
        
        for(j=0;j<EdofNum2;j++)//组装 C
            for(k=0;k<EdofNum2;k++)
                tripletList2.push_back(Eigen::Triplet<double>(dof2.TotalD[i][j],dof2.TotalD[i][k],LocC[j][k]));
        
        for(j=0;j<EdofNum1;j++)
        {
            delete [] LocA[j];delete [] LocAL2[j];delete [] LocAdiv[j];
        }
        for(j=0;j<EdofNum2;j++)    delete [] LocC[j];
        delete []LocA;delete []LocAL2; delete [] LocAdiv; delete [] LocC; delete []LocF;
    }
    
    StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
    StiffMatrixA.setFromTriplets(tripletList1.begin(),tripletList1.end());
    StiffMatrixC.setFromTriplets(tripletList2.begin(),tripletList2.end());
}

void DarcyStokesModel::GetSystem1()
{
    int i,j,k,EdofNum1,EdofNum2,N=0,N1=0,N2=0;
    for(i=0;i<ms.Element_Num;i++)
    {
        j=VE1.dof.Total_Num_PerElement[i]; k=VE2.dof.Total_Num_PerElement[i];
        N+=j*j*2; N1+=j*k; N2+=k*k;
    }
    std::vector< Eigen::Triplet<double> > tripletList; //for constructing sparsmatrix A
    tripletList.reserve(N);
    std::vector< Eigen::Triplet<double> > tripletList1; //for constructing sparsmatrix SMB
    tripletList1.reserve(N1);
    std::vector< Eigen::Triplet<double> > tripletList2; //for constructing sparsmatrix C
    tripletList2.reserve(N2);
    for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
    {
        EdofNum1=dof1.Total_Num_PerElement[i];EdofNum2=dof2.Total_Num_PerElement[i];
        double **LocA=new double*[EdofNum1];double **LocAL2=new double*[EdofNum1];
        double **LocAdiv=new double*[EdofNum1];
        double *LocF=new double[EdofNum1];double **LocC=new double*[EdofNum2];
        for(j=0;j<EdofNum1;j++)
        {
            LocA[j]=new double[EdofNum1];   LocAL2[j]=new double[EdofNum1];
            LocAdiv[j]=new double[EdofNum2];
        }
        for(j=0;j<EdofNum2;j++)    LocC[j]=new double[EdofNum2];
        MVE.GetMixedA(i,LocA);        MVE.GetMixedB(i,LocAdiv);
        VE1.GetA_L2(i, LocAL2);
        GetLocRHS(i,LocF);VE2.GetA_L2(i,LocC);
        for(j=0;j<EdofNum1;j++)//组装 A,B
        {
            RHS1[ dof1.TotalD[i][j] ]+=LocF[j];//组装右端项
            for(k=0;k<EdofNum1;k++)
                tripletList.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.TotalD[i][k],pde.epsilon*pde.epsilon*LocA[j][k]+LocAL2[j][k]));
            for(k=0;k<EdofNum2;k++)
                tripletList1.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof2.TotalD[i][k],LocAdiv[j][k]));
        }
        
        for(j=0;j<EdofNum2;j++)//组装 C
            for(k=0;k<EdofNum2;k++)
                tripletList2.push_back(Eigen::Triplet<double>(dof2.TotalD[i][j],dof2.TotalD[i][k],LocC[j][k]));
        
        for(j=0;j<EdofNum1;j++)
        {
            delete [] LocA[j];delete [] LocAL2[j];delete [] LocAdiv[j];
        }
        for(j=0;j<EdofNum2;j++)    delete [] LocC[j];
        delete []LocA;delete []LocAL2; delete [] LocAdiv; delete [] LocC; delete []LocF;
    }
    
    StiffMatrixA.setFromTriplets(tripletList.begin(),tripletList.end());
    SMB.setFromTriplets(tripletList1.begin(),tripletList1.end());
    StiffMatrixC.setFromTriplets(tripletList2.begin(),tripletList2.end());
}

double DarcyStokesModel::VelocityEnergyError()
{
    VectorXd X=uI.pVector-uh.pVector;
    double value=X.transpose()*StiffMatrixA*X;
    return sqrt(value);
}

double DarcyStokesModel::PressureL2Error()
{
    VectorXd X=pI.pVector-ph.pVector;
    double value=X.transpose()*StiffMatrixC*X;
    return sqrt(value);
}

double DarcyStokesModel::PressureMaxError()
{
    double value=0;
    int i;
    for(i=0;i<dof2.Dof_Num;i++)
        if(value<fabs(pI[i]-ph[i]))
            value=fabs(pI[i]-ph[i]);
    return value;
}

void DarcyStokesModel::OutputData()
{
    GetSystem();
    //deal with boundary condition for u
    int i,k,dofNum=dof1.Dof_Num+dof2.Dof_Num,BdofNum=ms.B_Node_Num*dof1.Num_PerNode+ms.B_Edge_Num*dof1.Num_PerEdge,
    IdofNum=dofNum-BdofNum;;
    int *Bdof=new int[BdofNum]; int *Idof=new int[IdofNum];int *IdofID=new int[dofNum];
    double *BdofVal=new double[BdofNum];
    for(i=0;i<BdofNum;i++)
        BdofVal[i]=0;
    VE1.GetBdof_BdofVal(pde.DirichletBoundary,Bdof,BdofVal);
    
    long N=StiffMatrix.nonZeros();
    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(N);
    
    VectorXd x0(dof1.Dof_Num+dof2.Dof_Num);
    x0.setZero();
    for(i=0;i<BdofNum;i++)
    {
        uh[Bdof[i]]=BdofVal[i]; x0(Bdof[i])=BdofVal[i];
    }
    RHS-=StiffMatrix*x0;
    
    FILE *fprow,*fpcol,*fpA,*fpB;
    fprow=fopen("./example/DSresult/Arow.txt","w");
    fpcol=fopen("./example/DSresult/Acol.txt","w");
    fpA=fopen("./example/DSresult/Aval.txt","w");
    fpB=fopen("./example/DSresult/B.txt","w");
    int col=0,row;
    for(i=0;i<dofNum;i++) //get the relation between dof and inner dof
        if(i<dof1.Dof_Num)
        {
            IdofID[i]=-1;
            if(dof1.DType[i]==0)
            {
                Idof[col]=i; IdofID[i]=col;  col++;
                fprintf(fpB,"%16.12f\n",RHS[i]);
            }
        }
        else
        {
            Idof[col]=i; IdofID[i]=col;  col++;
            fprintf(fpB,"%16.12f\n",RHS[i]);
        }
    
    for(col=0;col<IdofNum;col++)
    {
        i=Idof[col];
        for(SparseMatrix<double>::InnerIterator it(StiffMatrix,i);it;++it)
        {
            k=(int) it.row();
            if(k<dof1.Dof_Num)
            {
                if(dof1.DType[k]==0)
                {
                    row=IdofID[k];
                    tripletList.push_back(Eigen::Triplet<double>(row,col,it.value()));
                    fprintf(fprow,"%d ",row+1);fprintf(fpcol,"%d ",col+1);fprintf(fpA,"%16.12f ",it.value());
                }
            }
            else
            {
                row=IdofID[k];
                tripletList.push_back(Eigen::Triplet<double>(row,col,it.value()));
                fprintf(fprow,"%d ",row+1);fprintf(fpcol,"%d ",col+1);fprintf(fpA,"%16.12f ",it.value());
            }
        }
    }
    
    fclose(fprow); fclose(fpcol); fclose(fpA);fclose(fpB);
    fpB=fopen("./example/DSresult/X.txt","r");
    double temp;
    for(i=0;i<IdofNum;i++)
    {
        k=Idof[i];
        if(k<dof1.Dof_Num)
        {
            fscanf(fpB,"%lf",&temp);
            uh[k]=temp;
        }
        else
        {
            fscanf(fpB,"%lf",&temp);
            ph[k-dof1.Dof_Num]=temp;
        }
    }
    
    //add L_0^2 condition for ph
    double mean=ph.GetAverageValue();
    for(i=0;i<ph.ms.Element_Num;i++)
        ph[i*ph.dof.Num_PerElement]-=mean;
    
    delete []Bdof; delete []BdofVal; delete []Idof; delete []IdofID;
}
