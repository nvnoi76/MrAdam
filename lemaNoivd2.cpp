/*
LAMANOI by Nguyen Van Noi
Date : 11/01/2015
Modified : 13/07/2017
Version: 1.0.1
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 


#define HAVE_LAPACK
#define DBL_RAND_MAX (double)(RAND_MAX)
#define M_PI   3.14159265358979323846

//So tham so hoi quy
#define NUMPARAM  3

#define	_X0	p[0]
#define	_Y0	p[1]
#define	_R	p[2]


//So tham so bang du lieu
#define NUMCOL	  2

// So dong bang du lieu 
#define NUMLINE   82
#define NUMROW	  82

double gdata[NUMLINE][NUMCOL];

#define	_X		gdata[i][0]
#define	_Y		gdata[i][1]

double sqr(double x);
int ReadData(int n, char * Filename);
//Ham f(x)
void Tinhf(double* p, double* f, int numparam, int numrow)
{
	for (int i = 0; i < numrow; ++i)
		f[i] = sqr(_X-_X0)+sqr(_Y-_Y0) -sqr(_R);
}
void TinhJ(double* p, double* jac, int numparam, int numrow)
{
	int i, j;
	for (i = j = 0; i < numrow; ++i) 
	{
		jac[j++]=2*(_X0-_X);
		jac[j++]=2*(_Y0-_Y);
		jac[j++]=(-2)*_R;							
	}
}
void TinhG(double* p, double* g, int numparam, int numrow)
{
	double f[NUMROW];
	double jac[NUMROW][NUMPARAM];
	Tinhf(p,f,NUMPARAM,NUMROW);
	TinhJ(p,(double *)jac,NUMPARAM,NUMROW);
	int i,j;
	for(j=0;j<numparam;j++)
	{
		g[j]=0;
		for(i=0;i<numrow;i++)
			g[j]+=f[i]*jac[i][j];
	}
}
double TinhF(double* p, int numparam, int numrow)
{
	double f[NUMROW], F = 0;
	Tinhf(p,f,numparam,numrow);
	int i;
	for(i =0;i<numrow;i++)
	{
		F+= sqr(f[i]);
	}
	return F/2.0;
}
void xuatm(double g[], int m)
{
	for(int j=0;j<m;j++)
		printf("_a[%d]=  %.7g\n", j,g[j]);
	printf("\n");
}
#define ALPHA  0.01
#define BETA1  0.9 
#define BETA2  0.99

void Adam(double * p, int numparam, int numrow, int maxiter=300000,
		  double alpha = ALPHA, double beta1= BETA1, double beta2= BETA2, double eps=1e-8)
{
	double m[NUMPARAM],v[NUMPARAM], mhat[NUMPARAM], vhat[NUMPARAM], g[NUMPARAM];
	int j;
	for(j=0;j<numparam;j++)
		m[j] = v[j] = 0;
	double beta1t = beta1, beta2t = beta2, alphat = alpha, F=0;
    int i = 1;
	while (i<maxiter)
	{
		TinhG(p,g,numparam, numrow);
		for(j=0;j<numparam;j++)
		{
			m[j] = beta1t * m[j] + (1 - beta1t) * g[j];
            v[j] = beta2t * v[j] + (1 - beta2t) * g[j] * g[j];
            mhat[j] = m[j] / (1 - beta1t);
            vhat[j] = v[j] / (1 - beta2t);
		}
		beta1t = pow(beta1, i);
        beta2t = pow(beta2, i);
        alphat = alpha * sqrt(1 - beta2t) / (1 - beta1t);
		for(j=0;j<numparam;j++)
			p[j] = p[j] - alphat * mhat[j] / (sqrt(vhat[j])+ eps);
		i++;
		xuatm(g,numparam);

		F = TinhF(p,numparam,numrow);
		printf("i = %d => F =  %.7g\n",i, F);
	}
	xuatm(p,numparam);

}


void main()
{
	const int n= NUMROW, m= NUMPARAM;
	double p[NUMPARAM]={13, 15.5, 35};
	char * filein = "F:\\_BAIBAOLEMA_Adam\\AdamCircle\\datavd2.txt";	
	if (!ReadData(NUMLINE,filein))
		printf("Loi doc file %s!\n", filein), exit(-1);
	Adam(p,NUMPARAM, NUMROW);
	printf("X=[%.14g, %.14g, %.14g]T\n",_X0,_Y0,_R);
}
int ReadData(int n, char * Filename)
{
	register int i;	
	FILE * file = fopen(Filename,"rt");
	if (!file) 
		return 0;	
	for(i=0; i<n; ++i)
	{
		fscanf(file,"%lf",&_X);
		fscanf(file,"%lf",&_Y);
	}	
	fclose(file);
	
	file = fopen(".\\lemaNOIout1.csv","wt");
	if (!file)
		return 0;
	
	for(i=0; i<n; ++i)
	{
		fprintf(file," %.14g,",_X);
		fprintf(file," %.14g\n",_Y);
	}
	fclose(file);
	return 1;	
}
double sqr(double x)
{
	return x*x;
}