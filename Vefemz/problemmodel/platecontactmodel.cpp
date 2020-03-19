#include <iostream>
#include "problemmodel.h"


const double PlateContactProblem::nu=0.3;
void PlateContactProblem::Source(double x,double y,double *value)
{
    *value=3*(1-x*x)*(1-x*x)+3*(1-y*y)*(1-y*y)+(3*x*x-1)*(3*y*y-1);

}
void PlateContactProblem::FrictionalBoundary(double x,double y,double *value)
{

    if(x<1./6)
        *value=0.1;
    else
        *value=0.9;
 
//    *value=0.5;

}

void PlateContactProblem::ClampedBoundary(double x, double y, double *value)
{
    value[0]=0; value[1]=0;
}

void PlateContactProblem::Solution(double x, double y, double *value)
{
    value[0]=0;
}

void PlateContactModel::GetRHSL2B(int ElemID,double *LocFB)
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

void PlateContactModel::GetLocRHS(int ElemID,double *LocF)
{
    VE.GetRHS(ElemID,pde.Source,LocF);
}

void PlateContactModel::GetBInf()
{
    int i,en;
    double y,ey[2];
    double y1=ms.domain.y_min,y2=ms.domain.y_max;
    GaussLegendreQuadrature GLQ(VE.p+2);
    
    if(dof.Num_PerNode>0)
        for(i=0;i<ms.B_Node_Num;i++)
        {
            y=ms.Node[ ms.B_Node[i] ][1];
            if(fabs(y-y2)<0.000000001) //Gamma_D={y=1}
                BdofNum+=dof.Num_PerNode;
        }
    for(i=0;i<ms.B_Edge_Num;i++)
    {
        en=ms.Edge[ms.B_Edge[i] ][0];
        ey[0]=ms.Node[en][1];
        en=ms.Edge[ms.B_Edge[i] ][1];
        ey[1]=ms.Node[en][1];
        if(fabs(ey[0]-y1)<0.000000001&&fabs(ey[1]-y1)<0.000000001)
        {
            BQuadPtsNum+=GLQ.QuadPtsNum;
            ms.EdgeType[ms.B_Edge[i]]=2;
        }
    }
    if(dof.Num_PerEdge>0)
        for(i=0;i<ms.B_Edge_Num;i++)
        {
            en=ms.Edge[ms.B_Edge[i] ][0];
            ey[0]=ms.Node[en][1];
            en=ms.Edge[ms.B_Edge[i] ][1];
            ey[1]=ms.Node[en][1];
            if(fabs(ey[0]-y2)<0.000000001&&fabs(ey[1]-y2)<0.000000001)
                BdofNum+=dof.Num_PerEdge;
        }
}

void PlateContactModel::GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval)
{
    int i,j,k,m,n=0,en;
    double d,temp,value[2],x,y,ex[2],ey[2];
    double y2=ms.domain.y_max;
 
    if(dof.Num_PerNode>0)
        for(i=0;i<ms.B_Node_Num;i++)
        {
            x=ms.Node[ ms.B_Node[i] ][0];    y=ms.Node[ ms.B_Node[i] ][1];
            if(fabs(y-y2)<0.000000001) //Gamma_D={y=1}
            {
                BFunc(x,y,value);
                for(j=0;j<dof.Num_PerNode;j++)
                {
                    Bdof[n]=dof.NodeD[ ms.B_Node[i] ] [j];
                    Bdofval[n]=value[j];
                    n++;
                }
            }
            else
            {
                for(j=0;j<dof.Num_PerNode;j++)
                    dof.DType[ dof.NodeD[ ms.B_Node[i] ] [j] ]=2;
            }
        }
    
    if(dof.Num_PerEdge>0)
    {
        GaussLegendreQuadrature GLQ(VE.p+2);
        double jx,jy,stdx;
        for(i=0;i<ms.B_Edge_Num;i++)
        {
            en=ms.Edge[ms.B_Edge[i] ][0];
            ex[0]=ms.Node[en][0]; ey[0]=ms.Node[en][1];
            en=ms.Edge[ms.B_Edge[i] ][1];
            ex[1]=ms.Node[en][0]; ey[1]=ms.Node[en][1];
            if(fabs(ey[0]-y2)<0.000000001&&fabs(ey[1]-y2)<0.000000001)
            {
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
                        if(ElemType==0)
                        {
                            for(m=0;m<j%(dof.Num_PerEdge/2);m++)
                                temp*=stdx;
                            if(j<dof.Num_PerEdge/2)
                                Bdofval[n]+=temp*value[0]/2*GLQ.Weights[k];
                            else
                                Bdofval[n]+=temp*value[1]*d/2*GLQ.Weights[k];
                        }
                        else if(ElemType==1)
                        {
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
                    }
                    n++;
                }
            }
            else
            {
                for(j=0;j<dof.Num_PerEdge;j++)
                    dof.DType[ dof.EdgeD[ ms.B_Edge[i] ] [j] ]=2;
            }
        }
    }
    assert(n==BdofNum);
}

void PlateContactModel::GetSystem()
{
    int i,j,k,EdofNum,N=0;
    for(i=0;i<ms.Element_Num;i++)
    {
        j=VE.dof.Total_Num_PerElement[i];  N+=j*j;
    }
    std::vector< Eigen::Triplet<double> > tripletList; //for constructing sparsmatrix
    tripletList.reserve(N);
    
    for(i=0;i<ms.Element_Num;i++)//组装刚度矩阵和载荷
    {
        EdofNum=dof.Total_Num_PerElement[i];
        double **LocA=new double*[EdofNum];
        double *LocF=new double[EdofNum];
        for(j=0;j<EdofNum;j++)        LocA[j]=new double[EdofNum];
        VE.GetA_H2(i,LocA);
        GetLocRHS(i,LocF);
    
        for(j=0;j<EdofNum;j++)//组装
        {
            RHS[ dof.TotalD[i][j] ]+=LocF[j];//组装右端项
            for(k=0;k<EdofNum;k++)
                tripletList.push_back(Eigen::Triplet<double>(dof.TotalD[i][j],dof.TotalD[i][k],LocA[j][k]));
        }
        for(j=0;j<EdofNum;j++)
            delete [] LocA[j];
        delete []LocA; delete []LocF;
    }
    StiffMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
}

void PlateContactModel::InitialSolve(int *Idof,int *IdofID,VectorXd &b,VectorXd &x,SparseMatrix<double> &A)
{
    GetSystem();
    //deal with boundary condition
    int i,k,IdofNum=dof.Dof_Num-BdofNum;
    int *Bdof=new int[BdofNum];
    double *BdofVal=new double[BdofNum];
    for(i=0;i<BdofNum;i++)
        BdofVal[i]=0;
    GetBdof_BdofVal(pde.ClampedBoundary,Bdof,BdofVal);
 
    long N=StiffMatrix.nonZeros();
    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(N);
    
    for(i=0;i<BdofNum;i++)
        uh[Bdof[i]]=BdofVal[i];
    RHS-=StiffMatrix*uh.pVector;
    
    int col=0,row;
    for(i=0;i<dof.Dof_Num;i++) //get the relation between dof and inner dof
    {
        IdofID[i]=-1;
        if(dof.DType[i]!=1)
        {
            b[col]=RHS[i];
            Idof[col]=i;//inner dof to dof
            IdofID[i]=col;//dof to inner
            col++;
        }
    }
    
    for(col=0;col<IdofNum;col++)
    {
        i=Idof[col];
        for(SparseMatrix<double>::InnerIterator it(StiffMatrix,i);it;++it)
        {
            k=(int) it.row();
            if(dof.DType[k]!=1)
            {
                row=IdofID[k];
                tripletList.push_back(Eigen::Triplet<double>(row,col,it.value()));
            }
        }
    }

    A.setFromTriplets(tripletList.begin(),tripletList.end());
    //direct method
    Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    
    solver.compute(A);
    x=solver.solve(b);
    delete []Bdof; delete []BdofVal;
}

void PlateContactModel::IterativeSolve()
{
    double rho=.1,lam1,lam2,g[1];
    VectorXd lambda(BQuadPtsNum); lambda.setZero();
    int *vedof=new int[VE.p+1];double *ebasefun=new double[VE.p+1];
    int IdofNum=dof.Dof_Num-BdofNum;
    int *Idof=new int[IdofNum];int *IdofID=new int[dof.Dof_Num];
    VectorXd b(IdofNum),x(IdofNum);
    SparseMatrix<double> A(IdofNum,IdofNum);
    b.setZero();
    InitialSolve(Idof,IdofID,b,x,A);
    int i;
    GaussLegendreQuadrature GLQ(VE.p+2);
    Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    VectorXd R(dof.Dof_Num);
    solver.compute(A);
    
    for(int en=0;en<50;en++)
    {
        int gn=0;    lam1=0;lam2=0;
        for(i=0;i<IdofNum;i++)
            b[i]=RHS[Idof[i]];

        for(i=0;i<ms.B_Edge_Num;i++)
        {
            double x1=ms.Node[ ms.Edge[ ms.B_Edge[i] ][0] ][0];
            double x2=ms.Node[ ms.Edge[ ms.B_Edge[i] ][1] ][0];
            double y1=ms.Node[ ms.Edge[ ms.B_Edge[i] ][0] ][1];
            double y2=ms.Node[ ms.Edge[ ms.B_Edge[i] ][1] ][1];
            if(ms.EdgeType[ ms.B_Edge[i] ]==2)
            {
                double d=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
                vedof[0]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][0] ][0] ];
                vedof[1]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][1] ][0] ];
                vedof[2]=IdofID[ dof.EdgeD[ ms.B_Edge[i] ] [0] ];
                if(VE.p==3) vedof[3]=IdofID[ dof.EdgeD[ ms.B_Edge[i] ] [1] ];
                for(int k=0;k<GLQ.QuadPtsNum;k++)
                {
                    double stdx=GLQ.QuadPts[k]; //高斯积分点
                    double jx=(x1+x2)/2.+(x2-x1)*stdx/2., jy=(y1+y2)/2.+(y2-y1)*stdx/2.;
                    pde.FrictionalBoundary(jx, jy, g);
                    double temp=0;
                    if(VE.p==2)
                    {
                        ebasefun[0]=0.25*(-1-2*stdx+3*stdx*stdx);
                        ebasefun[1]=0.25*(-1+2*stdx+3*stdx*stdx);
                        ebasefun[2]=1.5*(1-stdx*stdx);
                        b[vedof[0]]-=g[0]*lambda[gn]*ebasefun[0]*d/2.0*GLQ.Weights[k];
                        b[vedof[1]]-=g[0]*lambda[gn]*ebasefun[1]*d/2.0*GLQ.Weights[k];
                        b[vedof[2]]-=g[0]*lambda[gn]*ebasefun[2]*d/2.0*GLQ.Weights[k];
                        temp=lambda[gn]+rho*g[0]*(x[vedof[0]]*ebasefun[0]+x[vedof[1]]*ebasefun[1]+x[vedof[2]]*ebasefun[2]);
                    }
                    else if(VE.p==3)
                    {
                        ebasefun[0]=0.25*(-1+3*stdx+3*stdx*stdx-5*stdx*stdx*stdx);
                        ebasefun[1]=0.25*(-1-3*stdx+3*stdx*stdx+5*stdx*stdx*stdx);
                        ebasefun[2]=1.5*(1-stdx*stdx);
                        ebasefun[3]=7.5*stdx*(1-stdx*stdx);
                        b[vedof[0]]-=g[0]*lambda[gn]*ebasefun[0]*d/2.0*GLQ.Weights[k];
                        b[vedof[1]]-=g[0]*lambda[gn]*ebasefun[1]*d/2.0*GLQ.Weights[k];
                        b[vedof[2]]-=g[0]*lambda[gn]*ebasefun[2]*d/2.0*GLQ.Weights[k];
                        b[vedof[3]]-=g[0]*lambda[gn]*ebasefun[2]*d/2.0*GLQ.Weights[k];
                        temp=lambda[gn]+rho*g[0]*(x[vedof[0]]*ebasefun[0]+x[vedof[1]]*ebasefun[1]+x[vedof[2]]*ebasefun[2]+x[vedof[3]]*ebasefun[3]);
                    }
                    else
                    {
                        assert(VE.p<4);
                    }
                    
                    if(temp>1) temp=1;
                    else if(temp<-1) temp=-1;
                    //        printf("en=%d   lam1=%f, lam2=%f\n",en,lambda[gn],rho*(mX[vedof[0]]*ebasefun[0]+mX[vedof[1]]*ebasefun[1]+mX[vedof[2]]*ebasefun[2]));
                    lam2+=fabs(temp-lambda[gn])*d/2.0*GLQ.Weights[k];
                    lam1+=fabs(lambda[gn])*d/2.0*GLQ.Weights[k];
                    lambda[gn]=temp;
                    gn++;
               
                }
            }
        }
        printf("en=%d   lam1=%f,lam2-lam1/lam1=%12.9f\n",en,lam1,lam2/lam1);
        if(lam2/lam1<0.00000001) break;
        x=solver.solve(b);
        
    }
    for(i=0;i<IdofNum;i++) uh[Idof[i]]=x[i];
    
    delete []vedof; delete []ebasefun; delete []Idof; delete []IdofID;
}

void PlateContactModel::IterativeSolve2()
{
    double rho=.1,lam1,lam2,g[1];
    double quadvalue[1];
    VectorXd lambda(BQuadPtsNum); lambda.setZero();
    int *vedof=new int[VE.p+1];double *ebasefun=new double[VE.p+1];
    int IdofNum=dof.Dof_Num-BdofNum;
    int *Idof=new int[IdofNum];
    int *IdofID=new int[dof.Dof_Num];
    VectorXd b(IdofNum),x(IdofNum);
    SparseMatrix<double> A(IdofNum,IdofNum);
    b.setZero();
    InitialSolve(Idof,IdofID,b,x,A);
    int i;
    GaussLegendreQuadrature GLQ(VE.p+2);
    Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    VectorXd R(dof.Dof_Num);
    solver.compute(A);

    for(int en=0;en<50;en++)
    {
        int gn=0;    lam1=0;lam2=0;
        for(i=0;i<IdofNum;i++)
            b[i]=RHS[Idof[i]];
        
        for(i=0;i<ms.B_Edge_Num;i++)
        {
            double x1=ms.Node[ ms.Edge[ ms.B_Edge[i] ][0] ][0];
            double x2=ms.Node[ ms.Edge[ ms.B_Edge[i] ][1] ][0];
            double y1=ms.Node[ ms.Edge[ ms.B_Edge[i] ][0] ][1];
            double y2=ms.Node[ ms.Edge[ ms.B_Edge[i] ][1] ][1];
            if(ms.EdgeType[ ms.B_Edge[i] ]==2)
            {
                double d=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
                vedof[0]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][0] ][0] ];
                vedof[1]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][1] ][0] ];
                if(VE.p==3) vedof[2]=IdofID[ dof.EdgeD[ ms.B_Edge[i] ] [0] ];
                uhEdgeQuadValue(ms.B_Edge[i],IdofID,x,quadvalue);
                for(int k=0;k<GLQ.QuadPtsNum;k++)
                {
                    double stdx=GLQ.QuadPts[k]; //高斯积分点
                    double jx=(x1+x2)/2.+(x2-x1)*stdx/2., jy=(y1+y2)/2.+(y2-y1)*stdx/2.;
                    pde.FrictionalBoundary(jx, jy, g);
                    double temp=0;
                    if(VE.p==2)
                    {
                        ebasefun[0]=0.25*(-1-2*stdx+3*stdx*stdx);
                        ebasefun[1]=0.25*(-1+2*stdx+3*stdx*stdx);
                        ebasefun[2]=1.5*(1-stdx*stdx);
                        b[vedof[0]]-=g[0]*lambda[gn]*ebasefun[0]*d/2.0*GLQ.Weights[k];
                        b[vedof[1]]-=g[0]*lambda[gn]*ebasefun[1]*d/2.0*GLQ.Weights[k];
                        double bvalue=g[0]*lambda[gn]*ebasefun[2]*d/2.0*GLQ.Weights[k];
                        AddtoRight(ms.B_Edge[i],bvalue,IdofID,b);
                        temp=lambda[gn]+rho*g[0]*(x[vedof[0]]*ebasefun[0]+x[vedof[1]]*ebasefun[1]+quadvalue[0]/d*ebasefun[2]);
                    }
                    else if(VE.p==3)
                    {
                        ebasefun[0]=0.25*(-1+3*stdx+3*stdx*stdx-5*stdx*stdx*stdx);
                        ebasefun[1]=0.25*(-1-3*stdx+3*stdx*stdx+5*stdx*stdx*stdx);
                        ebasefun[2]=1.5*(1-stdx*stdx);
                        ebasefun[3]=7.5*stdx*(1-stdx*stdx);
                        b[vedof[0]]-=g[0]*lambda[gn]*ebasefun[0]*d/2.0*GLQ.Weights[k];
                        b[vedof[1]]-=g[0]*lambda[gn]*ebasefun[1]*d/2.0*GLQ.Weights[k];
                        b[vedof[2]]-=g[0]*lambda[gn]*ebasefun[2]*d/2.0*GLQ.Weights[k];
                        double bvalue=g[0]*lambda[gn]*ebasefun[3]*d/2.0*GLQ.Weights[k];
                        AddtoRight(ms.B_Edge[i],bvalue,IdofID,b);
                        temp=lambda[gn]+rho*g[0]*(x[vedof[0]]*ebasefun[0]+x[vedof[1]]*ebasefun[1]+x[vedof[2]]*ebasefun[2]+quadvalue[0]/d*ebasefun[3]);
                    }
                    else
                    {
                        assert(VE.p<4);
                    }
                    
                    if(temp>1) temp=1;
                    else if(temp<-1) temp=-1;
                    //        printf("en=%d   lam1=%f, lam2=%f\n",en,lambda[gn],rho*(mX[vedof[0]]*ebasefun[0]+mX[vedof[1]]*ebasefun[1]+mX[vedof[2]]*ebasefun[2]));
                    lam2+=fabs(temp-lambda[gn])*d/2.0*GLQ.Weights[k];
                    lam1+=fabs(lambda[gn])*d/2.0*GLQ.Weights[k];
                    lambda[gn]=temp;
                    gn++;
                    
                }
            }
        }
        printf("en=%d   lam1=%f,lam2-lam1/lam1=%12.9f\n",en,lam1,lam2/lam1);
        if(lam2/lam1<0.00000001) break;
        x=solver.solve(b);
        
    }
    for(i=0;i<IdofNum;i++) uh[Idof[i]]=x[i];
  
    delete []vedof; delete []ebasefun; delete []Idof; delete []IdofID;
}

void PlateContactModel::IterativeSolve3()
{
    double rho=.1,lam1,lam2,g[1];
    double vertex1h,vertex2h,t[2];
    VectorXd lambda(BQuadPtsNum); lambda.setZero();
    int *vedof=new int[VE.p+4];double *ebasefun=new double[VE.p+2];
    int IdofNum=dof.Dof_Num-BdofNum;
    int *Idof=new int[IdofNum];int *IdofID=new int[dof.Dof_Num];
    VectorXd b(IdofNum),x(IdofNum);
    SparseMatrix<double> A(IdofNum,IdofNum);
    b.setZero();
    InitialSolve(Idof,IdofID,b,x,A);
    int i;
    GaussLegendreQuadrature GLQ(VE.p+2);
    Eigen::SimplicialLLT<SparseMatrix<double> > solver;
    VectorXd R(dof.Dof_Num);
    solver.compute(A);
    Getvcharlength();
 
    for(int en=0;en<50;en++)
    {
        int gn=0;    lam1=0;lam2=0;
        for(i=0;i<IdofNum;i++)
            b[i]=RHS[Idof[i]];
        
        for(i=0;i<ms.B_Edge_Num;i++)
        {
            double x1=ms.Node[ ms.Edge[ ms.B_Edge[i] ][0] ][0];
            double x2=ms.Node[ ms.Edge[ ms.B_Edge[i] ][1] ][0];
            double y1=ms.Node[ ms.Edge[ ms.B_Edge[i] ][0] ][1];
            double y2=ms.Node[ ms.Edge[ ms.B_Edge[i] ][1] ][1];
            if(ms.EdgeType[ ms.B_Edge[i] ]==2)
            {
                double d=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
                t[0]=(x2-x1)/d;t[1]=(y2-y1)/d;
                vertex1h=vcharlength[ ms.Edge[ ms.B_Edge[i] ][0] ];
                vertex2h=vcharlength[ ms.Edge[ ms.B_Edge[i] ][1] ];
                vedof[0]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][0] ][0] ];
                vedof[1]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][0] ][1] ];
                vedof[2]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][0] ][2] ];
                vedof[3]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][1] ][0] ];
                vedof[4]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][1] ][1] ];
                vedof[5]=IdofID[ dof.NodeD[ ms.Edge[ ms.B_Edge[i] ][1] ][2] ];
        
                for(int k=0;k<GLQ.QuadPtsNum;k++)
                {
                    double stdx=GLQ.QuadPts[k]; //高斯积分点
                    double jx=(x1+x2)/2.+(x2-x1)*stdx/2., jy=(y1+y2)/2.+(y2-y1)*stdx/2.;
                    pde.FrictionalBoundary(jx, jy, g);
                    double temp=0;
                    if(VE.p==2)
                    {
                        ebasefun[0]=0.25*(2-3*stdx+stdx*stdx*stdx);
                        ebasefun[1]=0.25*(1-stdx-stdx*stdx+stdx*stdx*stdx);
                        ebasefun[2]=0.25*(2+3*stdx-stdx*stdx*stdx);
                        ebasefun[3]=0.25*(-1-stdx+stdx*stdx+stdx*stdx*stdx);
                        b[vedof[0]]-=g[0]*lambda[gn]*ebasefun[0]*d/2.0*GLQ.Weights[k];
                        b[vedof[1]]-=g[0]*lambda[gn]*ebasefun[1]*d/2.0*GLQ.Weights[k]*t[0]/vertex1h;
                        b[vedof[2]]-=g[0]*lambda[gn]*ebasefun[1]*d/2.0*GLQ.Weights[k]*t[1]/vertex1h;
                        b[vedof[3]]-=g[0]*lambda[gn]*ebasefun[2]*d/2.0*GLQ.Weights[k];
                        b[vedof[4]]-=g[0]*lambda[gn]*ebasefun[3]*d/2.0*GLQ.Weights[k]*t[0]/vertex2h;
                        b[vedof[5]]-=g[0]*lambda[gn]*ebasefun[3]*d/2.0*GLQ.Weights[k]*t[1]/vertex2h;
                        temp=lambda[gn]+rho*g[0]*(x[vedof[0]]*ebasefun[0]+(x[vedof[1]]*t[0]+x[vedof[2]]*t[1])/vertex1h*ebasefun[1]
                        +x[vedof[3]]*ebasefun[2]+(x[vedof[4]]*t[0]+x[vedof[5]]*t[1])/vertex2h*ebasefun[3]);
                    }

                    if(temp>1) temp=1;
                    else if(temp<-1) temp=-1;
                    //        printf("en=%d   lam1=%f, lam2=%f\n",en,lambda[gn],rho*(mX[vedof[0]]*ebasefun[0]+mX[vedof[1]]*ebasefun[1]+mX[vedof[2]]*ebasefun[2]));
                    lam2+=fabs(temp-lambda[gn])*d/2.0*GLQ.Weights[k];
                    lam1+=fabs(lambda[gn])*d/2.0*GLQ.Weights[k];
                    lambda[gn]=temp;
               
                    gn++;
               
                }
            }
        }
   
        printf("en=%d   lam1=%f,lam2-lam1/lam1=%12.9f\n",en,lam1,lam2/lam1);
        if(lam2/lam1<0.00000001) break;
        x=solver.solve(b);
   
    }
    for(i=0;i<IdofNum;i++) uh[Idof[i]]=x[i];
  
    delete []vedof; delete []ebasefun; delete []Idof; delete []IdofID;
}

int PlateContactModel::FindElemID(int ElemID, PolyMesh &mesh)
{
    int i,j,k,EID=-1,nv=ms.ElementVertex_Num[ElemID],nv2,vflag;
    double x,y,x0,y0,x1,y1,x2,y2;
    for (i=0;i<mesh.Element_Num;i++)
    {
        nv2=mesh.ElementVertex_Num[i];
        x0=mesh.ElementBarycenter[i][0]; y0=mesh.ElementBarycenter[i][1];
        vflag=0;
        for(j=0;j<nv;j++)
        {
            x=ms.Node[ ms.Element[ElemID][j] ][0]; y=ms.Node[ ms.Element[ElemID][j] ][1];
            for(k=0;k<nv2;k++)
            {
                x1=mesh.Node[ mesh.Element[i][k] ][0]; y1=mesh.Node[ mesh.Element[i][k] ][1];
                x2=mesh.Node[ mesh.Element[i][(k+1)%nv2] ][0]; y2=mesh.Node[ mesh.Element[i][(k+1)%nv2] ][1];
                if((y1-y0)*(x-x0)-(x1-x0)*(y-y0)<=0&&(y2-y1)*(x-x1)-(x2-x1)*(y-y1)<=0&&(y0-y2)*(x-x2)-(x0-x2)*(y-y2)<=0)
                {
                    vflag++; break;
                }
            }
        }
        if(vflag==nv)
        {
            EID=i;break;
        }
    }
    
    return EID;
}

void PlateContactModel::uhEdgeQuadValue(int EdgeID,int *IdofID,VectorXd &xuh,double *value)
{
    int ElemID,dofID,i,j;
    if(ms.EdgeElement[EdgeID][0]!=-1)
        ElemID=ms.EdgeElement[EdgeID][0];
    else
        ElemID=ms.EdgeElement[EdgeID][1];
    int EdofNum=dof.Total_Num_PerElement[ElemID];
    double **Pistar=new double*[polydim],temp=0;
    for(i=0;i<polydim;i++)        Pistar[i]=new double[EdofNum];
    VE.GetPistar_H2(ElemID, Pistar);
    double x1=ms.Node[ ms.Edge[EdgeID][0] ][0];
    double x2=ms.Node[ ms.Edge[EdgeID][1] ][0];
    double y1=ms.Node[ ms.Edge[EdgeID][0] ][1];
    double y2=ms.Node[ ms.Edge[EdgeID][1] ][1];
    double d=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    GaussLegendreQuadrature GLQ(VE.p+2);
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(int k=0;k<GLQ.QuadPtsNum;k++)
    {
        double stdx=GLQ.QuadPts[k]; //高斯积分点
        double x=(x1+x2)/2.+(x2-x1)*stdx/2., y=(y1+y2)/2.+(y2-y1)*stdx/2.;
        for(i=0;i<EdofNum;i++)
        {
            dofID=dof.TotalD[ElemID][i];
            for(j=0;j<polydim;j++)
                temp+=xuh[ IdofID[dofID] ]*Pistar[j][i]*SMS.GetValue(j, x, y, hE, xE, yE)*pow(stdx,VE.p-2)*d/2.0*GLQ.Weights[k];
        }
    }
    *value=temp;
    for(i=0;i<polydim;i++)  delete []Pistar[i];
    delete []Pistar;

}
void PlateContactModel::AddtoRight(int EdgeID,double bvalue,int *IdofID,VectorXd &b)
{
    int ElemID,dofID,i,j;
    if(ms.EdgeElement[EdgeID][0]!=-1)
        ElemID=ms.EdgeElement[EdgeID][0];
    else
        ElemID=ms.EdgeElement[EdgeID][1];
    int EdofNum=dof.Total_Num_PerElement[ElemID];
    double **Pistar=new double*[polydim],temp=0;
    for(i=0;i<polydim;i++)        Pistar[i]=new double[EdofNum];
    VE.GetPistar_H2(ElemID, Pistar);
    double x1=ms.Node[ ms.Edge[EdgeID][0] ][0];
    double x2=ms.Node[ ms.Edge[EdgeID][1] ][0];
    double y1=ms.Node[ ms.Edge[EdgeID][0] ][1];
    double y2=ms.Node[ ms.Edge[EdgeID][1] ][1];
    GaussLegendreQuadrature GLQ(VE.p+2);
    double hE=ms.ElementDiameter[ElemID],xE=ms.ElementBarycenter[ElemID][0],yE=ms.ElementBarycenter[ElemID][1];
    for(i=0;i<EdofNum;i++)
    {
        temp=0;
        for(int k=0;k<GLQ.QuadPtsNum;k++)
        {
            double stdx=GLQ.QuadPts[k]; //高斯积分点
            double x=(x1+x2)/2.+(x2-x1)*stdx/2., y=(y1+y2)/2.+(y2-y1)*stdx/2.;
            for(j=0;j<polydim;j++)
                temp+=Pistar[j][i]*SMS.GetValue(j, x, y, hE, xE, yE)*pow(stdx,VE.p-2)/2.0*GLQ.Weights[k];
        }
        temp*=bvalue;
        dofID=dof.TotalD[ElemID][i];
        b[ IdofID[dofID] ]-=temp;        
    }
    
    for(i=0;i<polydim;i++)  delete []Pistar[i];
    delete []Pistar;
}

double PlateContactModel::EnergyError(PlateContactModel &pcmodel)
{
    double v1,v2,value=0;
    double x[3],y[3],hE,xE,yE,hE2,xE2,yE2;
    int i,j,k,m,n,nv,nv2,i2,EdofN,EdofN2,dofID;
    for(i=0;i<ms.Element_Num;i++)
    {
        i2=FindElemID(i, pcmodel.ms);
        nv=ms.ElementVertex_Num[i]; EdofN=dof.Total_Num_PerElement[i];
        nv2=pcmodel.ms.ElementVertex_Num[i2]; EdofN2=pcmodel.dof.Total_Num_PerElement[i2];
        double **PiH2=new double*[polydim];double **PiH22=new double*[polydim];
        for(j=0;j<polydim;j++)
        {
            PiH2[j]=new double[EdofN]; PiH22[j]=new double[EdofN2];
        }
        VE.GetPistar_H2(i, PiH2); pcmodel.VE.GetPistar_H2(i2, PiH22);
        hE2=pcmodel.ms.ElementDiameter[i2];xE2=pcmodel.ms.ElementBarycenter[i2][0];yE2=pcmodel.ms.ElementBarycenter[i2][1];
        hE=ms.ElementDiameter[i];xE=ms.ElementBarycenter[i][0];yE=ms.ElementBarycenter[i][1];
        double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
        for(j=0;j<nv;j++)
        {
            xpos[j]=ms.Node[ ms.Element[i][j] ][0]; ypos[j]=ms.Node[ ms.Element[i][j] ][1];
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
                v1=0;
                for(m=0;m<EdofN;m++)
                {
                    dofID=dof.TotalD[i][m];
                    for(n=0;n<polydim;n++)
                        v1+=uh[dofID]*PiH2[n][m]*SMS.GetDerivative(n,2,0, jx, jy, hE, xE, yE);
                }
                v2=0;
                for(m=0;m<EdofN2;m++)
                {
                    dofID=pcmodel.dof.TotalD[i2][m];
                    for(n=0;n<polydim;n++)
                        v2+=pcmodel.uh[dofID]*PiH22[n][m]*SMS.GetDerivative(n,2,0,jx, jy, hE2, xE2, yE2);
                }
                value+=(v1-v2)*(v1-v2)*vol_K*jw;//d^2/dx^2
                v1=0;
                for(m=0;m<EdofN;m++)
                {
                    dofID=dof.TotalD[i][m];
                    for(n=0;n<polydim;n++)
                        v1+=uh[dofID]*PiH2[n][m]*SMS.GetDerivative(n,1,1, jx, jy, hE, xE, yE);
                }
                v2=0;
                for(m=0;m<EdofN2;m++)
                {
                    dofID=pcmodel.dof.TotalD[i2][m];
                    for(n=0;n<polydim;n++)
                        v2+=pcmodel.uh[dofID]*PiH22[n][m]*SMS.GetDerivative(n,1,1,jx, jy, hE2, xE2, yE2);
                }
                value+=2*(v1-v2)*(v1-v2)*vol_K*jw;//d^2/dxy
                v1=0;
                for(m=0;m<EdofN;m++)
                {
                    dofID=dof.TotalD[i][m];
                    for(n=0;n<polydim;n++)
                        v1+=uh[dofID]*PiH2[n][m]*SMS.GetDerivative(n,0,2, jx, jy, hE, xE, yE);
                }
                v2=0;
                for(m=0;m<EdofN2;m++)
                {
                    dofID=pcmodel.dof.TotalD[i2][m];
                    for(n=0;n<polydim;n++)
                        v2+=pcmodel.uh[dofID]*PiH22[n][m]*SMS.GetDerivative(n,0,2,jx, jy, hE2, xE2, yE2);
                }
                value+=(v1-v2)*(v1-v2)*vol_K*jw;//d^2/dyy2
            }
        }
        
        for(j=0;j<polydim;j++)
        {
            delete []PiH2[j]; delete [] PiH22[j];
        }
        delete []PiH2; delete []PiH22; delete []xpos; delete []ypos;
    }
    
    return sqrt(value);
}

double PlateContactModel::EnergyErrorFB(PlateContactModel &pcmodel)
{
    double v1,v2,value=0,ymax=ms.domain.y_max;
    double x[3],y[3],hE,xE,yE,hE2,xE2,yE2;
    int i,j,k,m,n,nv,nv2,i2,EdofN,EdofN2,dofID,flag;
    for(i=0;i<ms.Element_Num;i++)
    {
        
        nv=ms.ElementVertex_Num[i];
        
        flag=0;
        for(j=0;j<nv;j++)
        {
            v1=ms.Node[ ms.Element[i][j] ][1];
            
            if(fabs(v1-ymax)<=0.000012) flag=1;
        }
        if(flag!=1) continue;
        
        EdofN=dof.Total_Num_PerElement[i];
        i2=FindElemID(i, pcmodel.ms);
        nv2=pcmodel.ms.ElementVertex_Num[i2]; EdofN2=pcmodel.dof.Total_Num_PerElement[i2];
        double **PiH2=new double*[polydim];double **PiH22=new double*[polydim];
        for(j=0;j<polydim;j++)
        {
            PiH2[j]=new double[EdofN]; PiH22[j]=new double[EdofN2];
        }
        VE.GetPistar_H2(i, PiH2); pcmodel.VE.GetPistar_H2(i2, PiH22);
        hE2=pcmodel.ms.ElementDiameter[i2];xE2=pcmodel.ms.ElementBarycenter[i2][0];yE2=pcmodel.ms.ElementBarycenter[i2][1];
        hE=ms.ElementDiameter[i];xE=ms.ElementBarycenter[i][0];yE=ms.ElementBarycenter[i][1];
        double *xpos=new double[nv];  double *ypos=new double[nv];//存储单元位置
        for(j=0;j<nv;j++)
        {
            xpos[j]=ms.Node[ ms.Element[i][j] ][0]; ypos[j]=ms.Node[ ms.Element[i][j] ][1];
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
                v1=0;
                for(m=0;m<EdofN;m++)
                {
                    dofID=dof.TotalD[i][m];
                    for(n=0;n<polydim;n++)
                        v1+=uh[dofID]*PiH2[n][m]*SMS.GetDerivative(n,2,0, jx, jy, hE, xE, yE);
                }
                v2=0;
                for(m=0;m<EdofN2;m++)
                {
                    dofID=pcmodel.dof.TotalD[i2][m];
                    for(n=0;n<polydim;n++)
                        v2+=pcmodel.uh[dofID]*PiH22[n][m]*SMS.GetDerivative(n,2,0,jx, jy, hE2, xE2, yE2);
                }
                value+=(v1-v2)*(v1-v2)*vol_K*jw;//d^2/dx^2
                v1=0;
                for(m=0;m<EdofN;m++)
                {
                    dofID=dof.TotalD[i][m];
                    for(n=0;n<polydim;n++)
                        v1+=uh[dofID]*PiH2[n][m]*SMS.GetDerivative(n,1,1, jx, jy, hE, xE, yE);
                }
                v2=0;
                for(m=0;m<EdofN2;m++)
                {
                    dofID=pcmodel.dof.TotalD[i2][m];
                    for(n=0;n<polydim;n++)
                        v2+=pcmodel.uh[dofID]*PiH22[n][m]*SMS.GetDerivative(n,1,1,jx, jy, hE2, xE2, yE2);
                }
                value+=2*(v1-v2)*(v1-v2)*vol_K*jw;//d^2/dxy
                v1=0;
                for(m=0;m<EdofN;m++)
                {
                    dofID=dof.TotalD[i][m];
                    for(n=0;n<polydim;n++)
                        v1+=uh[dofID]*PiH2[n][m]*SMS.GetDerivative(n,0,2, jx, jy, hE, xE, yE);
                }
                v2=0;
                for(m=0;m<EdofN2;m++)
                {
                    dofID=pcmodel.dof.TotalD[i2][m];
                    for(n=0;n<polydim;n++)
                        v2+=pcmodel.uh[dofID]*PiH22[n][m]*SMS.GetDerivative(n,0,2,jx, jy, hE2, xE2, yE2);
                }
                value+=(v1-v2)*(v1-v2)*vol_K*jw;//d^2/dyy2
            }
        }
        
        for(j=0;j<polydim;j++)
        {
            delete []PiH2[j]; delete [] PiH22[j];
        }
        delete []PiH2; delete []PiH22; delete []xpos; delete []ypos;
    }
    
    return sqrt(value);
}

void PlateContactModel::Getvcharlength()
{
    int i,j,vID;
    vcharlength=new double[ms.Node_Num];
    int *venum=new int[ms.Node_Num];
    for(i=0;i<ms.Node_Num;i++)
    {
        vcharlength[i]=0;venum[i]=0;
    }
    for(i=0;i<ms.Element_Num;i++)
        for(j=0;j<ms.ElementVertex_Num[i];j++)
        {
            vID=ms.Element[i][j];
            vcharlength[vID]+=ms.ElementDiameter[i];
            venum[vID]++;
        }
    for(i=0;i<ms.Node_Num;i++)
    {
        vcharlength[i]/=venum[i];
    }
    delete []venum;
}

