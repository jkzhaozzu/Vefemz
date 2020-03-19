//
//  polynomialspace.cpp
//  VEM2DforMac
//
//  Created by 张蓓 on 2019/1/3.
//  Copyright © 2019年 Jikun Zhao. All rights reserved.
//
#include <stdio.h>
#include <iostream>
#include "polynomialspace.h"

double ScaledMonomialSpace::GetValue(int i, double x, double y, double hE, double xE, double yE)
{
    value[0]=1.;
    for(int j=1;j<=p;j++)
    {
        int lastj=(j-1)*j/2;
        int jstart=j*(j+1)/2;
        for(int k=0;k<j;k++)
            value[jstart+k]=value[lastj+k]*(x-xE)/hE;
        value[jstart+j]=value[lastj+j-1]*(y-yE)/hE;
    }
    return value[i];
}


double ScaledMonomialSpace::GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE)
{
	if(p-(dxm+dyn)<0)	return 0;
    GetValue(i,x,y,hE,xE,yE);
	double vi;
	int i_row=0,i_col=0,last,j;
	for(j=0;j<=p;j++)	
		if(i<(j+2)*(j+1)/2)
		{
			i_row=j;	i_col=i-(j+1)*j/2;
			break;
		}

	if(i_row-i_col<dxm||i_col<dyn)	return 0;
	else
	{
		last=(i_row-1-(dxm+dyn)+2)*(i_row-1-(dxm+dyn)+1)/2;
		vi=value[last+(i_col-dyn)];
		for(j=i_row-i_col;j>i_row-i_col-dxm;j--)
			vi*=j/hE;
		for(j=i_col;j>i_col-dyn;j--)
			vi*=j/hE;
		return vi;
	}
}

void ScaledMonomialSpaceV::GetValue(int i, double x, double y, double hE, double xE, double yE,double *value)
{
	int m=i/2;
	if(i%2==0)
	{
		value[0]=SMS.GetValue(m,x,y,hE,xE,yE);
		value[1]=0;
	}
	else
	{
		value[0]=0;
		value[1]=SMS.GetValue(m,x,y,hE,xE,yE);
	}
}

void ScaledMonomialSpaceV::GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE,double *value)
{
	int m=i/2;
	if(i%2==0)
	{
		value[0]=SMS.GetDerivative(m,dxm,dyn,x,y,hE,xE,yE);
		value[1]=0;
	}
	else
	{
		value[0]=0;
		value[1]=SMS.GetDerivative(m,dxm,dyn,x,y,hE,xE,yE);
	}
}

void GradPolynomialSpace::GetValue(int i, double x, double y, double hE, double xE, double yE,double *value)
{
    value[0]=SMS.GetDerivative(i+1,1,0,x,y,hE,xE,yE);
	value[1]=SMS.GetDerivative(i+1,0,1,x,y,hE,xE,yE);	
}

void GradPolynomialSpace::GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE,double *value)
{
	value[0]=SMS.GetDerivative(i+1,dxm+1,dyn,x,y,hE,xE,yE);
	value[1]=SMS.GetDerivative(i+1,dxm,dyn+1,x,y,hE,xE,yE);
}

void GradPolynomialOrthogonalSpace::GetValue(int i, double x, double y, double hE, double xE, double yE,double *value)
{
    value[0]=SMS.GetValue(i,x,y,hE,xE,yE)*(y-yE)/hE;
	value[1]=-SMS.GetValue(i,x,y,hE,xE,yE)*(x-xE)/hE;	
}

void GradPolynomialOrthogonalSpace::GetDerivative(int i,int dxm,int dyn,double x,double y,double hE,double xE,double yE,double *value)
{
	if(dyn==0)
		value[0]=SMS.GetDerivative(i,dxm,dyn,x,y,hE,xE,yE)*(y-yE)/hE;
	else
		value[0]=SMS.GetDerivative(i,dxm,dyn,x,y,hE,xE,yE)*(y-yE)/hE+dyn*SMS.GetDerivative(i,dxm,dyn-1,x,y,hE,xE,yE)/hE;
	if(dxm==0)
		value[1]=-SMS.GetDerivative(i,dxm,dyn,x,y,hE,xE,yE)*(x-xE)/hE;	
	else
		value[1]=-SMS.GetDerivative(i,dxm,dyn,x,y,hE,xE,yE)*(x-xE)/hE-dxm*SMS.GetDerivative(i,dxm-1,dyn,x,y,hE,xE,yE)/hE;	
}
