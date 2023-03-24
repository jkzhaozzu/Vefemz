#include <iostream>
using namespace std;
#include "stdlib.h"
#include "stdio.h"
#include <iomanip>
#include "time.h"

#include "mesh.h"
#include "dof.h"
#include "polynomialspace.h"
#include "quadrature.h"
#include "mixedvirtualelement.h"
#include "problemmodel.h"



void Curl4VEHcurl2()
{
    FILE* fpdata=fopen("./example/Curl4Result/error.txt","w");
    fprintf(fpdata,"h         ||e||        ||curle||       ||curl^2e||\n\n");
    for(int meshID=1;meshID<4;meshID++)
    {
        
        clock_t start,end;
        double t1=0,t2=0;
        int k=2;//k>=2
        FILE *fp;
      //  fp=GetMeshPointer(trimesh,meshID);
      //  fp=GetMeshPointer(rectmesh,meshID);
        fp=GetMeshPointer(unipolymesh,meshID);
        PolyMesh ms(fp);
        fclose(fp);
        
        VEHcurl2 VE1(k,ms);
        VEPkC0 VE2(k-1,ms);
        MVEHcurl2 MVE(VE1,VE2);
        Curl4Model curl4VEM(MVE);
        
        start=clock();
        curl4VEM.Solve();
 
        end=clock();
        double uhEnergy=curl4VEM.uI.Norm(curl4VEM.Curl2A),uhCurl=curl4VEM.uI.Norm(curl4VEM.CurlA);
        double uhL2=curl4VEM.uI.Norm(curl4VEM.AL2);
        double uEnergyError=curl4VEM.uEnergyError(),uCurlError=curl4VEM.uCurlError(),uL2Error=curl4VEM.uL2Error();
      
        
        t1=double(end-start)/CLOCKS_PER_SEC;
        cout<< "meshID= "<<meshID<<endl;
        cout<<"solving time: " << t1 <<" s"<< endl;
        cout<<"uhEnergy="<<uhEnergy<<" uEnergyErr=  "<<uEnergyError<<" relative err "<<uEnergyError/uhEnergy<<endl;
        cout<<"uhCurl="<<uhCurl<<" uCurlErr=  "<<uCurlError<<" relative err "<<uCurlError/uhCurl<<endl;
        cout<<"uhL2="<<uhL2<<" uL2Err=  "<<uL2Error<<" relative err "<<uL2Error/uhL2<<endl;
        double psiError=curl4VEM.psiH1Error();
        cout<<"psiH1Erro="<<psiError<<endl;
        fprintf(fpdata,"%f  %E  %E  %E\n",ms.ElementDiameter[0],uL2Error,uCurlError,uEnergyError);
        
        //      cout<<stokes.SMA<<endl;
        
        //        cout<<"pI.aver="<<stokes.pI.GetAverageValue()<<" ph.aver="<<stokes.ph.GetAverageValue()<<endl;
        //        for(int i=0;i<stokes.dof2.Dof_Num;i++)
        //            cout<<"i= "<<i<<"  "<<stokes.pI[i]<<"  "<<stokes.ph[i]<<"  "<<endl;

        
  /*      int ElemID=0;
         int i,j,m,NdofE=dof1.Total_Num_PerElement[ElemID],NdofE2=dof2.Total_Num_PerElement[ElemID], polyVdim=(k+2)*(k+1),pdim=(k+1)*k/2,polydim=k*(k-1);
        
         double ** GH1=new double*[pdim];double ** GL2=new double*[polydim],**Gtest=new double*[pdim];
         double ** BDH1=new double*[pdim];double ** BL2=new double*[polydim];
         double ** BH1=new double*[pdim];double ** BDL2=new double*[polydim];
         double ** D=new double*[NdofE];double ** D2=new double*[NdofE2];double ** MixedB=new double*[NdofE];
        double ** AH1=new double*[NdofE];double ** G2=new double*[NdofE];double ** BDG2=new double*[NdofE];
         for(i=0;i<pdim;i++)
         {
         GH1[i]=new double[pdim];Gtest[i]=new double[polyVdim];
         BH1[i]=new double[NdofE];BDH1[i]=new double[polyVdim];
         }
        for(i=0;i<polydim;i++)
        {GL2[i]=new double[polyVdim];BL2[i]=new double[NdofE];BDL2[i]=new double[NdofE];}
         for(i=0;i<NdofE;i++)
         {
             D[i]=new double[polyVdim];
             AH1[i]=new double[NdofE];MixedB[i]=new double[NdofE2];G2[i]=new double[pdim];BDG2[i]=new double[pdim];
         }
        for(i=0;i<NdofE2;i++) D2[i]=new double[pdim];
   
      //   VE1.GetPistar_L2(ElemID,BH1);
        MVE.GetMixedB(ElemID, MixedB);

         for(i=0;i<NdofE;i++)
         {
             for(j=0;j<NdofE2;j++)
             {      
                         cout<<MixedB[i][j]<<"  ";
    //             BDL2[i][j]=0;
    //             for(m=0;m<NdofE2;m++)
  //                   BDL2[i][j]+=BL2[i][m]*D2[m][j];
  //               cout<<"("<<i<<","<<j<<"):"<<G2[i][j]-BDL2[i][j]<<endl;
                 
             }
             cout<<endl;
         }
 */

        
    }
    fclose(fpdata);
}


