#include <iostream>
#include "problemmodel.h"
#include "mathfunction.h"

void NavierStokesProblem::Source1(double x,double y,double *value)
{
    double v[6];
    Velocity(x,y,v);
    value[0]=-8*PI*PI*cos(2*PI*x)*sin(2*PI*y)+4*PI*PI*sin(2*PI*y)+y*y+v[0]*v[1]+v[3]*v[2];
    value[1]=8*PI*PI*cos(2*PI*y)*sin(2*PI*x)-4*PI*PI*sin(2*PI*x)+2*x*y+v[0]*v[4]+v[3]*v[5];
    //    value[0]=-8*PI*PI*cos(2*PI*x)*sin(2*PI*y)+4*PI*PI*sin(2*PI*y)+cos(x)+v[0]*v[1]+v[3]*v[2];
    //    value[1]=8*PI*PI*cos(2*PI*y)*sin(2*PI*x)-4*PI*PI*sin(2*PI*x)-cos(y)+v[0]*v[4]+v[3]*v[5];
    
}

void NavierStokesProblem::Source2(double x,double y,double *value)
{
    value[0]=0;
}

void NavierStokesProblem::DirichletBoundary(double x,double y,double *value)
{
    value[0]=0; value[1]=0;
}
void NavierStokesProblem::Velocity(double x,double y,double *value)
{
    value[0]=-cos(2*PI*x)*sin(2*PI*y)+sin(2*PI*y);
    value[1]=2*PI*sin(2*PI*x)*sin(2*PI*y);
    value[2]=2*PI*cos(2*PI*y)*(1-cos(2*PI*x));
    value[3]=sin(2*PI*x)*cos(2*PI*y)-sin(2*PI*x);
    value[4]=2*PI*cos(2*PI*x)*(cos(2*PI*y)-1);
    value[5]=-2*PI*sin(2*PI*x)*sin(2*PI*y);
    /*
     value[0]=-cos(2*PI*x)*sin(2*PI*y)+sin(2*PI*y);
     value[1]=2*PI*sin(2*PI*x)*sin(2*PI*y);
     value[2]=2*PI*cos(2*PI*y)*(1-cos(2*PI*x));
     value[3]=sin(2*PI*x)*cos(2*PI*y)-sin(2*PI*x);
     value[4]=2*PI*cos(2*PI*x)*(cos(2*PI*y)-1);
     value[5]=-2*PI*sin(2*PI*x)*sin(2*PI*y);
     */
}

void NavierStokesProblem::Pressure(double x,double y,double *value)
{
    value[0]=x*y*y-1./6;    value[1]=0;    value[2]=0;
    //    value[0]=sin(x)-sin(y);    value[1]=0;    value[2]=0;
}

void NavierStokesModel::GetRHSL2B(int ElemID,double *LocFB)
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

void NavierStokesModel::GetLocRHS(int ElemID,double *LocF)
{
    VE1.GetRHS(ElemID,pde.Source1,LocF);
    /*    int i,j,EdofNum=dof1.Total_Num_PerElement[ElemID];
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
     */
}

void NavierStokesModel::GetMatrixRHS()
{
    ;
    int i,j,k,EdofNum1,EdofNum2,N=0,N1=0,N2=0;
    for(i=0;i<ms.Element_Num;i++)
    {
        j=VE1.dof.Total_Num_PerElement[i]; k=VE2.dof.Total_Num_PerElement[i];
        N+=j*k; N1+=j*j; N2+=k*k;
    }
    std::vector< Eigen::Triplet<double> > tripletListAdiv;
    tripletListAdiv.reserve(N);
    std::vector< Eigen::Triplet<double> > tripletListA;
    tripletListA.reserve(N1);
    std::vector< Eigen::Triplet<double> > tripletListAL2;
    tripletListAL2.reserve(N1);
    std::vector< Eigen::Triplet<double> > tripletListC;
    tripletListC.reserve(N2);
    for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
    {
        EdofNum1=dof1.Total_Num_PerElement[i];EdofNum2=dof2.Total_Num_PerElement[i];
        double **LocA=new double*[EdofNum1];double **LocAL2=new double*[EdofNum1];double **LocAdiv=new double*[EdofNum1];
        double **LocC=new double*[EdofNum2];double *LocF=new double[EdofNum1];        for(j=0;j<EdofNum1;j++)
        {
            LocA[j]=new double[EdofNum1]; LocAL2[j]=new double[EdofNum1]; LocAdiv[j]=new double[EdofNum2];
        }
        for(j=0;j<EdofNum2;j++)    LocC[j]=new double[EdofNum2];
        MVE.GetMixedA(i,LocA); VE1.GetA_L2(i,LocAL2); VE2.GetA_L2(i,LocC); MVE.GetMixedB(i,LocAdiv);GetLocRHS(i,LocF);
        for(j=0;j<EdofNum1;j++)//组装 A,B
        {
            RHS[ dof1.TotalD[i][j] ]+=LocF[j];//组装右端项
            for(k=0;k<EdofNum1;k++)
            {
                tripletListA.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.TotalD[i][k],LocA[j][k]));
                tripletListAL2.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.TotalD[i][k],LocAL2[j][k]));
            }
            for(k=0;k<EdofNum2;k++)
            {
                tripletListAdiv.push_back(Eigen::Triplet<double>(dof2.TotalD[i][k],dof1.TotalD[i][j],LocAdiv[j][k]));
            }
        }
        for(j=0;j<EdofNum2;j++)//组装 C
            for(k=0;k<EdofNum2;k++)
                tripletListC.push_back(Eigen::Triplet<double>(dof2.TotalD[i][j],dof2.TotalD[i][k],LocC[j][k]));
        
        for(j=0;j<EdofNum1;j++)
        {
            delete [] LocA[j];delete [] LocAL2[j]; delete [] LocAdiv[j];
        }
        
        for(j=0;j<EdofNum2;j++)    delete [] LocC[j];
        delete []LocA; delete []LocAL2; delete [] LocAdiv; delete [] LocC; delete [] LocF;
    }
    
    MatrixA.setFromTriplets(tripletListA.begin(),tripletListA.end());
    MatrixAL2.setFromTriplets(tripletListAL2.begin(),tripletListAL2.end());
    MatrixAdivT.setFromTriplets(tripletListAdiv.begin(),tripletListAdiv.end());
    MatrixC.setFromTriplets(tripletListC.begin(),tripletListC.end());
    
}
void NavierStokesModel::GetProjectVal(int BasisID,double x,double y,double hE,double xE,double yE, double **A,double *value)
{
    int i,polydim=(VE1.p+2)*(VE1.p+1);
    double basisvalue[2];
    value[0]=0;value[1]=0;
    for(i=0;i<polydim;i++)
    {
        SMSV.GetValue(i,x,y,hE,xE,yE,basisvalue);
        value[0]+=A[i][BasisID]*basisvalue[0];value[1]+=A[i][BasisID]*basisvalue[1];
    }
}

void NavierStokesModel::GetProjectGradVal(int BasisID,double x,double y,double hE,double xE,double yE, double **A,double *value)
{
    int i,polydim=(VE1.p+2)*(VE1.p+1);
    double basisvaluedx[2],basisvaluedy[2];
    value[0]=0;value[1]=0;value[2]=0;value[3]=0;
    for(i=0;i<polydim;i++)
    {
        SMSV.GetDerivative(i,1,0,x,y,hE,xE,yE,basisvaluedx);
        SMSV.GetDerivative(i,0,1,x,y,hE,xE,yE,basisvaluedy);
        value[0]+=A[i][BasisID]*basisvaluedx[0];value[1]+=A[i][BasisID]*basisvaluedy[0];
        value[2]+=A[i][BasisID]*basisvaluedx[1];value[3]+=A[i][BasisID]*basisvaluedy[1];
    }
    
}
void NavierStokesModel::GetProjectGradVal(int BasisID,double x,double y,double hE,double xE,double yE,double **A1,double **A2,double *value)
{
    int i,polydim=(VE1.p+1)*VE1.p;
    double basisvaluedx[2],basisvaluedy[2];
    value[0]=0;value[1]=0;value[2]=0;value[3]=0;
    for(i=0;i<polydim;i++)
    {
        SMSV.GetValue(i,x,y,hE,xE,yE,basisvaluedx);
        SMSV.GetValue(i,x,y,hE,xE,yE,basisvaluedy);
        value[0]+=A1[i][BasisID]*basisvaluedx[0];value[1]+=A2[i][BasisID]*basisvaluedy[0];
        value[2]+=A1[i][BasisID]*basisvaluedx[1];value[3]+=A2[i][BasisID]*basisvaluedy[1];
    }
}

void NavierStokesModel::GetLocNonLinearMat(int ElemID,VEMFunction &uh0,double **LocA)
{
    int j,k,k1,k2,m,nv,EdofNum1,polydim=(VE1.p+2)*(VE1.p+1),pdim=(VE1.p+1)*VE1.p;
    double x[3],y[3],vigrad[4],vj[2],vuh[2];
    EdofNum1=dof1.Total_Num_PerElement[ElemID];nv=ms.ElementVertex_Num[ElemID];
    double *locuh0=new double[EdofNum1];
    double **LocPiH1=new double*[polydim];double **LocPiL2=new double*[polydim];
    double **LocdxPiL2=new double*[pdim]; double **LocdyPiL2=new double*[pdim];
    for(j=0;j<polydim;j++)
    {
        LocPiH1[j]=new double[EdofNum1];    LocPiL2[j]=new double[EdofNum1];
        if(j<pdim)
        {
            LocdxPiL2[j]=new double[EdofNum1]; LocdyPiL2[j]=new double[EdofNum1];
        }
    }
    for(j=0;j<EdofNum1;j++)
    {
        locuh0[j]=uh0[dof1.TotalD[ElemID][j]];
        for(k=0;k<EdofNum1;k++)
            LocA[j][k]=0;
    }
    
    VE1.GetPistar_H1(ElemID,LocPiH1);     VE1.GetPistar_L2(ElemID,LocPiL2);
    VE1.GetdxPistar_L2(ElemID,LocdxPiL2);VE1.GetdyPistar_L2(ElemID,LocdyPiL2);
    
    double *xpos=new double[nv]; double *ypos=new double[nv];//存储单元位置
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(j=0;j<nv;j++)
    {
        xpos[j]=ms.Node[ ms.Element[ElemID][j] ][0];   ypos[j]=ms.Node[ ms.Element[ElemID][j] ][1];
    }
    for(j=0;j<nv;j++)
    {
        x[0]=xE;x[1]=xpos[j];x[2]=xpos[(j+1)%nv];
        y[0]=yE;y[1]=ypos[j];y[2]=ypos[(j+1)%nv];
        double vol_K=0.5*fabs((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]));
        for(k=0;k<TQ.QuadPtsNum;k++)
        {
            double j1=TQ.QuadPts[k][0]; double j2=TQ.QuadPts[k][1]; double j3=TQ.QuadPts[k][2];     //积分点
            double jx=x[0]*j1+x[1]*j2+x[2]*j3; double jy=y[0]*j1+y[1]*j2+y[2]*j3;
            double jw=TQ.Weights[k];//积分权重
            for(k1=0;k1<EdofNum1;k1++)
                for(k2=0;k2<EdofNum1;k2++)
                {
                    //GetProjectGradVal(k1,jx,jy,hE,xE,yE,LocPiH1,vigrad);
                    GetProjectGradVal(k1,jx,jy,hE,xE,yE,LocdxPiL2,LocdyPiL2,vigrad);
                    GetProjectVal(k2,jx,jy,hE,xE,yE,LocPiL2,vj);
                    double sum[2];
                    sum[0]=0;sum[1]=0;
                    for(m=0;m<EdofNum1;m++)
                    {
                        GetProjectVal(m,jx,jy,hE,xE,yE,LocPiL2,vuh);
                        sum[0]+=locuh0[m]*vuh[0];sum[1]+=locuh0[m]*vuh[1];
                    }
                    LocA[k1][k2]+=0.5*((sum[0]*vigrad[0]+sum[1]*vigrad[1])*vj[0]+(sum[0]*vigrad[2]+sum[1]*vigrad[3])*vj[1])*vol_K*jw;
                    
                    GetProjectVal(k1,jx,jy,hE,xE,yE,LocPiL2,vj);
                    //    GetProjectGradVal(k2,jx,jy,hE,xE,yE,LocPiH1,vigrad);
                    GetProjectGradVal(k2,jx,jy,hE,xE,yE,LocdxPiL2,LocdyPiL2,vigrad);
                    LocA[k1][k2]+=0.5*((sum[0]*vigrad[0]+sum[1]*vigrad[1])*vj[0]+(sum[0]*vigrad[2]+sum[1]*vigrad[3])*vj[1])*vol_K*jw;
                }
        }
    }
    
    for(j=0;j<polydim;j++)
    {
        delete [] LocPiH1[j];delete [] LocPiL2[j];
        if(j<pdim)
        {
            delete []LocdxPiL2[j]; delete LocdyPiL2[j];
        }
        
    }
    delete []LocPiH1; delete []LocPiL2; delete []xpos; delete []ypos;delete []locuh0;    delete []LocdxPiL2; delete []LocdyPiL2;
}
void NavierStokesModel::Solve()
{
    int i,j,k,mi,EdofNum1,EdofNum2,dofNum1=dof1.Dof_Num,dofNum2=dof2.Dof_Num;
    GetMatrixRHS();
    long N=MatrixA.nonZeros()+2*MatrixAdivT.nonZeros();
    for(i=0;i<ms.Element_Num;i++)
    {
        j=VE1.dof.Total_Num_PerElement[i]; N+=j*j;
    }
    for(mi=0;mi<20;mi++)
    {
        std::vector< Eigen::Triplet<double> > tripletList;
        tripletList.reserve(N);
        for(j=0;j<dofNum1;j++)
        {
            for(SparseMatrix<double>::InnerIterator it(MatrixA,j);it;++it)
                tripletList.push_back(Eigen::Triplet<double>(it.row(),j,it.value()));
            for(SparseMatrix<double>::InnerIterator it(MatrixAdivT,j);it;++it)
            {
                tripletList.push_back(Eigen::Triplet<double>(it.row()+dofNum1,j,it.value()));
                tripletList.push_back(Eigen::Triplet<double>(j,it.row()+dofNum1,it.value()));
            }
        }
        for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
        {
            EdofNum1=dof1.Total_Num_PerElement[i];
            double **LocA=new double*[EdofNum1];
            for(j=0;j<EdofNum1;j++)
                LocA[j]=new double[EdofNum1];
            GetLocNonLinearMat(i,uh0,LocA);
            for(j=0;j<EdofNum1;j++)//组装
            {
                
                for(k=0;k<EdofNum1;k++)
                    tripletList.push_back(Eigen::Triplet<double>(dof1.TotalD[i][j],dof1.TotalD[i][k],LocA[j][k]));
            }
            for(j=0;j<EdofNum1;j++)
                delete []LocA[j];
            delete [] LocA;
        }
        StiffMatrix.setZero();
        StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
        SolveSystem();
        if(IterationError()<0.000001) break;
        
        for(i=0;i<dofNum1;i++)    uh0[i]=uh[i];
        for(i=0;i<dofNum2;i++)    ph0[i]=ph[i];
    }
}

void NavierStokesModel::SolveSystem()
{
    int i,j,k,dofNum=dof1.Dof_Num+dof2.Dof_Num,BdofNum=ms.B_Node_Num*dof1.Num_PerNode+ms.B_Edge_Num*dof1.Num_PerEdge,
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
    x0=RHS-StiffMatrix*x0;
    
    VectorXd b(IdofNum),x(IdofNum);
    int col=0,row;
    for(i=0;i<dofNum;i++) //get the relation between dof and inner dof
        if(i<dof1.Dof_Num)
        {
            IdofID[i]=-1;
            if(dof1.DType[i]==0)
            {
                b[col]=x0[i];  Idof[col]=i; IdofID[i]=col;  col++;
            }
        }
        else
        {
            b[col]=x0[i];  Idof[col]=i; IdofID[i]=col;  col++;
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
    
    // iterative method
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    //    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
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

double NavierStokesModel::IterationError()
{
    double value1,value2;
    VectorXd X=uh0.pVector-uh.pVector;
    value1=X.transpose()*MatrixAL2*X;
    VectorXd Y=ph0.pVector-ph.pVector;
    value2=Y.transpose()*MatrixC*Y;
    
    std::cout<<"iteration error="<<(sqrt(value1)+sqrt(value2))<<std::endl;
    return (sqrt(value1)+sqrt(value2));
}

double NavierStokesModel::VelocityEnergyError()
{
    double value;
    VectorXd X=uI.pVector-uh.pVector;
    value=X.transpose()*MatrixA*X;
    return sqrt(value);
}

double NavierStokesModel::PressureL2Error()
{
    double value;
    VectorXd X=pI.pVector-ph.pVector;
    value=X.transpose()*MatrixC*X;
    return sqrt(value);
}

double NavierStokesModel::PressureMaxError()
{
    double value=0;
    int i;
    for(i=0;i<dof2.Dof_Num;i++)
        if(value<fabs(pI[i]-ph[i]))
            value=fabs(pI[i]-ph[i]);
    return value;
}

