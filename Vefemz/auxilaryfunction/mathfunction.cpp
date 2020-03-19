#include <iostream>

int factorial(int m)
{
	if(m==0) return 1;
	int v=1;
	for(int i=1;i<=m;i++)
		v*=i;
	return v;
}

//Gauss消元法,主对角占优
void GaussSolve(int n,double **AA,double *bb,double *X)
{
    double** A=new double*[n];
	double* b=new double[n];
	int i,j,k,l;
	for(k=0;k<n;k++)
		A[k]=new double[n];
	for(k=0;k<n;k++)
	{
		b[k]=bb[k];
		for(l=0;l<n;l++)
			A[k][l]=AA[k][l];
	}
	double t;
	for(k=0;k<n-1;k++)
		for(i=k+1;i<n;i++)
		{
			t=-A[i][k]/A[k][k];
			for(int j=k+1;j<n;j++)
				A[i][j]+=t*A[k][j];
			b[i]+=t*b[k];
		}
	
	for(i=n-1;i>=0;i--)
	{
        
		X[i]=b[i]/A[i][i];
		for(j=i+1;j<n;j++)
			X[i]-=A[i][j]*X[j]/A[i][i];
        
	}
	for(k=0;k<n;k++)
		delete[] A[k];
	delete[] A;
	delete[] b;
	
}

int Sign(double x)
{
	if(x>0)
		return 1;
	else if(x==0)
		return 0;
	else 
		return -1;
}

