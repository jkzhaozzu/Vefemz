#ifndef poissonmodel_h
#define poissonmodel_h

#include "virtualelement.h"
#include "mixedvirtualelement.h"
#include "vemfunction.h"

/*
PoissonProblem is
-Delta u=f in Omega
u=g on boundary
*/
class PoissonProblem
{
public:
	PoissonProblem()
	{ 
		sourcedim=1; valuedim=3;	
	}
	~PoissonProblem(){;};
	static void Source(double x,double y,double *value);	//f	
	static void DirichletBoundary(double x,double y,double *value); //g
	static void Solution(double x,double y,double *value); //u

	int sourcedim;// dim of source
	int valuedim; // dim of value of u and grad u
	
};

// Poisson solver
class PoissonModel
{
public:
	PoissonModel(VirtualElement &ve):SMS(ve.SMS),TQ(ve.TQ),ms(ve.ms),dof(ve.dof),VE(ve),uh(ve),uI(ve),StiffMatrix(dof.Dof_Num,dof.Dof_Num),RHS(dof.Dof_Num)
	{
        polydim=(VE.p+2)*(VE.p+1)/2; RHS.setZero();
		uI.Interpolate(pde.Solution);		
	}
	~PoissonModel(){};
	void GetLocRHS(int ElemID,double *LocF);//get local RHS on element ElemID
	void Solve();//deal with BD and get uh
	void GetSystem();// get StiffMatrix and RHS
	double EnergyError();
	double MaxError();

	int polydim;
	PolyMesh &ms;
    DegreeofFreedom &dof;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
	VirtualElement &VE;
	VEMFunction uh;//numerical solution
	VEMFunction uI;//interpolation of solution
	PoissonProblem pde;
    SparseMatrix<double> StiffMatrix;
    VectorXd RHS;

};


/*ParabolicProblem is
du/dt-Delta u=f in Omega, t\in (0,T)
u=g on boundary,  t\in (0,T)
u(x,0)=u_0 in Omega
*/

class ParabolicProblem
{
public:
	ParabolicProblem()
	{ 
		sourcedim=1; valuedim=3;	
	}
	~ParabolicProblem(){;};
	static void Source(double x,double y,double t,double *value);	//f	
	static void DirichletBoundary(double x,double y,double t,double *value); //g
	static void InitialData(double x,double y,double *value); //u0  
	static void Solution(double x,double y,double t,double *value); //u

	int sourcedim;// dim of source
	int valuedim; // dim of value of u and grad u
};

// ParabolicProblem solver
class ParabolicModel
{
public:
	ParabolicModel(VirtualElement &ve,int dtnum):SMS(ve.SMS),TQ(ve.TQ),ms(ve.ms),dof(ve.dof),VE(ve),uh(ve),uh0(ve),uI(ve),
    StiffMatrix(dof.Dof_Num,dof.Dof_Num),MatrixL2(dof.Dof_Num,dof.Dof_Num),MatrixH1(dof.Dof_Num,dof.Dof_Num),RHS(dof.Dof_Num)
	{
		polydim=(VE.p+2)*(VE.p+1)/2;	T=1; dt=T/dtnum;dtn=dtnum;
        tn=1; RHS.setZero();
		uh0.Interpolate(pde.InitialData);
		uI.Interpolate(pde.Solution,tn);

	}

	~ParabolicModel(){};
	void GetLocRHS(int ElemID,double t,double *LocF);//get local RHS on element ElemID at time t
	void GetMatrix();
	void GetLV(double t);// get load vector (RHS) at time t
	void Solve();// backward diference method for time
	double ErrorH1();
	double ErrorL2();

	int polydim,dtn;
	double tn0,tn,dt,T;
	PolyMesh &ms;
    DegreeofFreedom &dof;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
	VirtualElement &VE;
	VEMFunction uh;//numerical solution at tn
	VEMFunction uh0;//numerical solution at last step
	VEMFunction uI;//interpolation of solution at t=tn
	ParabolicProblem pde;
    SparseMatrix<double> StiffMatrix;
    SparseMatrix<double> MatrixL2;
    SparseMatrix<double> MatrixH1;
    VectorXd RHS;

};


/*
Elasticity Problem is
 \sigma=C\epsilon(u) in Omega
-div\sigma(u)=f in Omega
u=g on boundary
*/

class ElasticityProblem
{
public:
	ElasticityProblem()
	{ 
		sourcedim=2; valuedim=6;
	}
	~ElasticityProblem(){;};
	static void Source(double x,double y,double *value);	//f	
	static void DirichletBoundary(double x,double y,double *value); //g
	static void Solution(double x,double y,double *value); //u
    static void Stress(double x,double y,double *value); //sigma

	int sourcedim;// dim of source
	int valuedim; // dim of value of u and grad u
	static const double mu;
	static double const lambda;
};

// Elasticity solver
class ElasticityModel
{
public:
	ElasticityModel(VirtualElement &ve):SMS(ve.SMS),TQ(ve.TQ),ms(ve.ms),dof(ve.dof),VE(ve),uh(ve),uI(ve),StiffMatrix(dof.Dof_Num,dof.Dof_Num),RHS(dof.Dof_Num)
	{
        polydim=(VE.p+2)*(VE.p+1); RHS.setZero();
		uI.Interpolate(pde.Solution);		
	}

	~ElasticityModel(){};
    void GetRHSL2B(int ElemID,double *LocFB);//get the quadrature of RHS with space SMS on element ElemID
    void GetLocRHS(int ElemID,double *LocF);//get local RHS on element ElemID
    void Solve();//deal with BD and get uh
    void GetSystem();// get StiffMatrix and RHS
	double EnergyError();
	double MaxError();

	int polydim;
	PolyMesh &ms;
    DegreeofFreedom &dof;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
	VirtualElement &VE;
	VEMFunction uh;//numerical solution
	VEMFunction uI;//interpolation of solution
	ElasticityProblem pde;
    SparseMatrix<double> StiffMatrix;
    VectorXd RHS;
};


/*
Stokes Problem is
-\Delta u+grad p=f, div u=0 in Omega
u=u_0 on boundary
*/


class StokesProblem
{
public:
	StokesProblem()
	{ 
		sourcedim=2; valuedim=6;
	}
	~StokesProblem(){;};
	static void Source1(double x,double y,double *value);	//f	
	static void Source2(double x,double y,double *value);	//g
	static void DirichletBoundary(double x,double y,double *value); //u
	static void Velocity(double x,double y,double *value); //u
	static void Pressure(double x,double y,double *value); //p

	int sourcedim;// dim of source f
	int valuedim; // dim of value of u and grad u
};

//Stokes solver
class StokesModel
{
public:
    StokesModel(MixedVirtualElement &mve):
    MVE(mve),VE1(mve.VE1),VE2(mve.VE2),SMSV(VE1.SMS),SMS(VE2.SMS),TQ(VE1.TQ),ms(VE1.ms),dof1(VE1.dof),dof2(VE2.dof),
    uh(VE1),uI(VE1),ph(VE2),pI(VE2),StiffMatrix(dof1.Dof_Num+dof2.Dof_Num,dof1.Dof_Num+dof2.Dof_Num),StiffMatrixA(dof1.Dof_Num,dof1.Dof_Num),
    StiffMatrixC(dof2.Dof_Num,dof2.Dof_Num),RHS(dof1.Dof_Num+dof2.Dof_Num)
    {
        polydim=(VE1.p+2)*(VE1.p+1); RHS.setZero();
        uI.Interpolate(pde.Velocity);    pI.Interpolate(pde.Pressure);
    };
 
	~StokesModel(){};
	void GetRHSL2B(int ElemID,double *LocFB);
	void GetLocRHS(int ElemID,double *LocF);//get local RHS for f on element ElemID
	void GetSystem();
	void Solve();
	double VelocityEnergyError();
	double PressureL2Error();
	double PressureMaxError();

	int polydim;   
	MixedVirtualElement &MVE;
	VirtualElement &VE1;
	VirtualElement &VE2;
	PolynomialSpace &SMSV;
	PolynomialSpace &SMS;
	TriangleQuadrature &TQ;
	PolyMesh &ms;
    DegreeofFreedom &dof1;
	DegreeofFreedom &dof2;
	VEMFunction uh;//numerical solution
	VEMFunction uI;//interpolation of solution
	VEMFunction ph;//numerical solution
	VEMFunction pI;//interpolation of solution
	StokesProblem pde;
	SparseMatrix<double> StiffMatrix,StiffMatrixA,StiffMatrixC;
	VectorXd RHS;
};

/*
 Darcy-Stokes Problem is
 u-\epsilon*\Delta u+grad p=f, in Omega
 div u=0 in Omega
 u=u_0 on boundary
 */


class DarcyStokesProblem
{
public:
    DarcyStokesProblem()
    {
        sourcedim=2; valuedim=6;
    }
    ~DarcyStokesProblem(){;};
    static void Source1(double x,double y,double *value);    //f
    static void Source2(double x,double y,double *value);    //g
    static void DirichletBoundary(double x,double y,double *value); //u
    static void Velocity(double x,double y,double *value); //u
    static void Pressure(double x,double y,double *value); //p
    
    int sourcedim;// dim of source f
    int valuedim; // dim of value of u and grad u
    static const double epsilon;
};

//Darcy-Stokes solver
class DarcyStokesModel
{
public:
    DarcyStokesModel(MixedVirtualElement &mve):
    MVE(mve),VE1(mve.VE1),VE2(mve.VE2),SMSV(VE1.SMS),SMS(VE2.SMS),TQ(VE1.TQ),ms(VE1.ms),dof1(VE1.dof),dof2(VE2.dof),
    uh(VE1),uI(VE1),ph(VE2),pI(VE2),StiffMatrix(dof1.Dof_Num+dof2.Dof_Num,dof1.Dof_Num+dof2.Dof_Num),StiffMatrixA(dof1.Dof_Num,dof1.Dof_Num),
    StiffMatrixC(dof2.Dof_Num,dof2.Dof_Num),RHS(dof1.Dof_Num+dof2.Dof_Num),SMB(dof1.Dof_Num,dof2.Dof_Num),RHS1(dof1.Dof_Num)
    {
        polydim=(VE1.p+2)*(VE1.p+1); RHS.setZero();
        uI.Interpolate(pde.Velocity);    pI.Interpolate(pde.Pressure);
    };
    
    ~DarcyStokesModel(){};
    void GetRHSL2B(int ElemID,double *LocFB);
    void GetLocRHS(int ElemID,double *LocF);//get local RHS for f on element ElemID
    void GetSystem();
    void GetSystem1();
    void Solve();
    void UzawaSolve();
    void OutputData();
    double VelocityEnergyError();
    double PressureL2Error();
    double PressureMaxError();
    
 
    int polydim;
    MixedVirtualElement &MVE;
    VirtualElement &VE1;
    VirtualElement &VE2;
    PolynomialSpace &SMSV;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
    PolyMesh &ms;
    DegreeofFreedom &dof1;
    DegreeofFreedom &dof2;
    VEMFunction uh;//numerical solution
    VEMFunction uI;//interpolation of solution
    VEMFunction ph;//numerical solution
    VEMFunction pI;//interpolation of solution
    DarcyStokesProblem pde;
    SparseMatrix<double> StiffMatrix,StiffMatrixA,SMB,StiffMatrixC;
    VectorXd RHS,RHS1;
};


/*
Navier-Stokes Problem is
-\Delta u+(u\cdot\grad) u+grad p=f, div u=0 in Omega
u=u_0 on boundary
*/


class NavierStokesProblem
{
public:
	NavierStokesProblem()
	{ 
		sourcedim=2; valuedim=6;
	}
	~NavierStokesProblem(){;};
	static void Source1(double x,double y,double *value);	//f	
	static void Source2(double x,double y,double *value);	//g
	static void DirichletBoundary(double x,double y,double *value); //u
	static void Velocity(double x,double y,double *value); //u
	static void Pressure(double x,double y,double *value); //p

	int sourcedim;// dim of source f
	int valuedim; // dim of value of u and grad u
};

//Navier-Stokes solver
class NavierStokesModel
{
public:
    
    NavierStokesModel(MixedVirtualElement &mve):MVE(mve),VE1(mve.VE1),VE2(mve.VE2),SMSV(VE1.SMS),SMS(VE2.SMS),TQ(VE1.TQ),ms(VE1.ms),dof1(VE1.dof),dof2(VE2.dof),
    uh0(VE1),uh(VE1),uI(VE1),ph0(VE2),ph(VE2),pI(VE2),StiffMatrix(dof1.Dof_Num+dof2.Dof_Num,dof1.Dof_Num+dof2.Dof_Num),
    MatrixA(dof1.Dof_Num,dof1.Dof_Num),MatrixAL2(dof1.Dof_Num,dof1.Dof_Num),MatrixAdivT(dof2.Dof_Num,dof1.Dof_Num),
    MatrixC(dof2.Dof_Num,dof2.Dof_Num),RHS(dof1.Dof_Num+dof2.Dof_Num)
    {
        polydim=(VE1.p+2)*(VE1.p+1); RHS.setZero();
        uI.Interpolate(pde.Velocity);    pI.Interpolate(pde.Pressure);
    };

	~NavierStokesModel(){};
	void GetRHSL2B(int ElemID,double *LocFB);
	void GetLocRHS(int ElemID,double *LocF);//get local RHS for f on element ElemID
	void GetProjectVal(int BasisID,double x,double y,double hE,double xE,double yE,double **A,double *value);
	void GetProjectGradVal(int BasisID,double x,double y,double hE,double xE,double yE,double **A,double *value);
	void GetProjectGradVal(int BasisID,double x,double y,double hE,double xE,double yE,double **A1,double **A2,double *value);
	void GetLocNonLinearMat(int ElemID,VEMFunction &uh0,double **LocA);
	void GetMatrixRHS();//get StiffMatrixA,StiffMatrixAL2,StiffMatrixC,RHS
	void Solve();
    void SolveSystem();
	double IterationError();
	double VelocityEnergyError();
//	double VelocityL2Error();
	double PressureL2Error();
	double PressureMaxError();
	
	int polydim;   
	MixedVirtualElement &MVE;
	VirtualElement &VE1;
	VirtualElement &VE2;
	PolynomialSpace &SMSV;
	PolynomialSpace &SMS;
	TriangleQuadrature &TQ;
	PolyMesh &ms;
    DegreeofFreedom &dof1;
	DegreeofFreedom &dof2;
	VEMFunction uh0;//initial numerical solution
	VEMFunction uh;//numerical solution
	VEMFunction uI;//interpolation of solution
	VEMFunction ph0;//initial numerical solution
	VEMFunction ph;//numerical solution
	VEMFunction pI;//interpolation of solution
	StokesProblem pde;
    SparseMatrix<double> StiffMatrix,MatrixA,MatrixAL2,MatrixAdivT,MatrixC;//H1 matrix MatrixA for uh,L2 MatrixC for ph
    VectorXd RHS;
};


/*
BiharmonicProblem is
Delta^2 u=f in Omega
u=g on boundary
du/dn=h on boundary
*/
class BiharmonicProblem
{
public:
	BiharmonicProblem()
	{ 
		sourcedim=1; valuedim=3;	
	}
	~BiharmonicProblem(){;};
	static void Source(double x,double y,double *value);	//f	
	static void DirichletBoundary(double x,double y,double *value); //g
    static void DirichletBoundary2(double x,double y,double *value); //g
	static void Solution(double x,double y,double *value); //u

	int sourcedim;// dim of source
	int valuedim; // dim of value of u and grad u
	
};

// Biharmonic solver
class BiharmonicModel
{
public:
	BiharmonicModel(VirtualElement &ve):SMS(ve.SMS),TQ(ve.TQ),ms(ve.ms),dof(ve.dof),VE(ve),uh(ve),uI(ve),StiffMatrix(dof.Dof_Num,dof.Dof_Num),RHS(dof.Dof_Num)
	{
        polydim=(VE.p+2)*(VE.p+1)/2; RHS.setZero();
		uI.Interpolate(pde.Solution);		
	}

	~BiharmonicModel(){};
	void GetLocRHS(int ElemID,double *LocF);//get local RHS on element ElemID
	void GetRHSL2B(int ElemID,double *LocFB);//get the quadrature of RHS with space P_k-2 on element ElemID
    void Solve();//deal with BD and get uh
    void Solve2();//deal with BD and get uh
    void GetSystem();// get StiffMatrix and RHS
	double EnergyError();
	double MaxError();

	int polydim;
	PolyMesh &ms;
    DegreeofFreedom &dof;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
	VirtualElement &VE;
	VEMFunction uh;//numerical solution
	VEMFunction uI;//interpolation of solution
	BiharmonicProblem pde;
    SparseMatrix<double> StiffMatrix;
    VectorXd RHS;
};

/*
 Plate contact problem is
 a(u,v-u)+j(v)-j(u) >= (f,v-u), \forall v\in V=H^1_\Gamma_D(\Omega)
 where a(u,v)=(\Delta u,\Delta v)+(1-nu)(2(u12,v12)-(u11,v22)-(u22,v11),
 j(v)=<g,|v|>_\Gamma_C
 0<\nu<0.5 0<g<1
 \Gamma_C: the frictional contact boundary
 \Gamma_D: the clamped boundary
 \Gamma_N: the free boundary
 */

class PlateContactProblem
{
public:
    PlateContactProblem()
    {
        sourcedim=1; valuedim=3;
    }
    ~PlateContactProblem(){;};
    static void Source(double x,double y,double *value);    //f
    static void FrictionalBoundary(double x,double y,double *value); //g on Gamma_C
    static void ClampedBoundary(double x,double y,double *value); //on Gamma_D
    static void Solution(double x,double y,double *value); //u
    
    int sourcedim;// dim of source
    int valuedim; // dim of value of u and grad u
    static const double nu;
};

//Plate contact solver
class PlateContactModel
{
public:
    PlateContactModel(VirtualElement &ve):SMS(ve.SMS),TQ(ve.TQ),ms(ve.ms),dof(ve.dof),VE(ve),uh(ve),uI(ve),StiffMatrix(dof.Dof_Num,dof.Dof_Num),RHS(dof.Dof_Num)
    {
        polydim=(VE.p+2)*(VE.p+1)/2; RHS.setZero();BdofNum=0; BQuadPtsNum=0;
        ElemType=0;
        GetBInf();
    }
    PlateContactModel(VirtualElement &ve,int elemtype):SMS(ve.SMS),TQ(ve.TQ),ms(ve.ms),dof(ve.dof),VE(ve),uh(ve),uI(ve),StiffMatrix(dof.Dof_Num,dof.Dof_Num),RHS(dof.Dof_Num)
    {
        polydim=(VE.p+2)*(VE.p+1)/2; RHS.setZero();BdofNum=0; BQuadPtsNum=0;
        ElemType=elemtype;
        GetBInf();
    }
    
    ~PlateContactModel(){delete []vcharlength; };
    void GetLocRHS(int ElemID,double *LocF);//get local RHS on element ElemID
    void GetRHSL2B(int ElemID,double *LocFB);//get the quadrature of RHS with space P_k-2 on element ElemID
    void GetBdof_BdofVal(FunctionP BFunc,int *Bdof,double *Bdofval);// get boundary dof ID and corresponding dof value on Gamma_d
    void GetBInf();//get the number BdofNum of dof on boundary \Gamma_D and BQuadPtsNum
    void InitialSolve(int *Idof,int *IdofID,VectorXd &b,VectorXd &x,SparseMatrix<double> &A);//deal with BD and get uh0
    void IterativeSolve();//deal with BD and get uh for VEC0H2NC
    void IterativeSolve2();//deal with BD and get uh for VEMorleyType
    void IterativeSolve3();//deal with BD and get uh for VEPkH2
    void GetSystem();// get StiffMatrix and RHS
    int FindElemID(int ElemID,PolyMesh &mesh); //find mesh element ID in which the element ElemID
    
    //for VEMorleyType
    void uhEdgeQuadValue(int EdgeID,int *IdofID,VectorXd &xuh,double *value);// get some moments up to order k-2
    void AddtoRight(int EdgeID,double bvalue,int *IdofID,VectorXd &b);

    double EnergyError(PlateContactModel &pcmodel);//compute error w.r.t. numerical solution on coarse mesh
    double EnergyErrorFB(PlateContactModel &pcmodel);//compute error w.r.t. numerical solution on coarse mesh near frictional boundary
//    double MaxError();

    int polydim,BdofNum,BQuadPtsNum;//BQuadPtsNum: quadrature points number on boundary Gamma_C
    PolyMesh &ms;
    DegreeofFreedom &dof;
    PolynomialSpace &SMS;
    TriangleQuadrature &TQ;
    VirtualElement &VE;
    int ElemType;//initial 0 for VEC0H2NC,1 for VEMoreyType
    VEMFunction uh;//numerical solution
    VEMFunction uI;//interpolation of solution
    PlateContactProblem pde;
    SparseMatrix<double> StiffMatrix;
    VectorXd RHS;
private:
    void Getvcharlength();
    double *vcharlength;//the element diameter average w.r.t. each vertex
    
};
#endif
