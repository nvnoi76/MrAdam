#include "MrAdam.h"
double gdata[NUMROW][NUMCOL];
double teta = 0.01;
double sqr(double x)
{
	return x*x;
}
int ReadData(const char * Filename, int numrow, int numcol, double data[NUMROW][NUMCOL])
{

	FILE* file = fopen(Filename, "rt");
	if (!file) 
		return 0;
	//MR	PA	WOPT	WW	XM3	SS	QU	PI	LL	N200	XMD  MAU    SUBMAU   ID
	int i, j;
	for (i = 0; i < numrow; i++)
		for(j=0; j< numcol; j++)
			fscanf(file, "%lf", &data[i][j]);
	fclose(file);
	file = fopen("dataout.csv", "wt");
	if (!file)		return 0;
	for (i = 0; i < numrow; ++i)
	{
		for(j=0; j< numcol; j++)
			fprintf(file, " %.14g,", data[i][j]);
		fprintf(file, "\n", ID);
	}
	fclose(file);
	return 1;
}
// Tinh he so k1
double Tinhk1(int i, double* p)
{
	double k1;
	k1 = (A11 + A12 * (WOPT - WW) / WOPT) * pow(XM3, A2) + A3 * pow(SS / 100.0, A4) + A5 * QU + A6 * PI + A7 * (LL - WW) + A8 * (WOPT - WW) + A9 * (N200 - A10);
	return k1;
}
//Tinh co so cua luy thua k2
double Tinhcoso(int i, double* p)
{
	double coso;
	coso = 9.0 * PA * ((1.0 / (3.0 * XMD)) + (XM3 / (XMD * XMD))) / 2.0;
	return coso;
}
// Tinh he so k2
double Tinhk2(int i, double* p)
{
	double k2;
	k2 = (B11 + B12 * (WW - WOPT)) * pow(XM3, B2) + B3 * pow((SS / 100.0), B4) + B5 * pow(QU, B6) + B7 * PI + B8 * LL;
	return k2;
}
// Tinh gia tri f = Mr/Pa = k1* (coso ^ k2);
void Tinhf(double* p, double* f, int numparam, int numrow)
{
	int i;
	double k1, k2, coso;
	for (i = 0; i < numrow; ++i)
	{
		k1 = Tinhk1(i, p);
		coso = Tinhcoso(i, p);
		k2 = Tinhk2(i, p);
		f[i] = k1 * pow(coso, k2);
	}
}
/*
  Tinh ma tran Jacobian cho 20 tham so a11,a12,a2 ... 10, b11, b12, b2, ... b8
  theo cong thuc tinh dao ham rieng        df/da = df/dk * dk/da
  f = k1 * coso ^ k2 = k1 * mu
  df/da = df/dk1 * dk1/da = coso ^ k2 * dk1/da
						  = mu * dk1/da
  df/db = df/dk2 * dk2/db = k1 * coso ^ k2 * lncoso * dk2/db
						  = k1* mu * lncoso * dk2/db
						  = x * dk2/db
*/
void TinhJ(double* p, double* jac, int numparam, int numrow)
{
	int i, j;
	double k1, k2, coso, x, mu;
	for (i = j = 0; i < numrow; ++i) {
		k1 = Tinhk1(i, p);
		coso = Tinhcoso(i, p);
		k2 = Tinhk2(i, p);
		mu = pow(coso, k2);
		x = k1 * mu * log(coso);
		//tinh dao ham rieng cho a11
		jac[j++] = mu * pow(XM3, A2);
		//tinh dao ham rieng cho a12
		jac[j++] = mu * (WOPT - WW) * pow(XM3, A2) / WOPT;
		//tinh dao ham rieng cho a2
		jac[j++] = mu * (A11 + A12 * (WOPT - WW) / WOPT) * pow(XM3, A2) * log(XM3);
		//tinh dao ham rieng cho a3
		jac[j++] = mu * pow(SS / 100.0, A4);
		//tinh dao ham rieng cho a4
		jac[j++] = mu * A3 * pow(SS / 100.0, A4) * log(SS / 100.0);
		//tinh dao ham rieng cho a5
		jac[j++] = mu * QU;
		//tinh dao ham rieng cho a6
		jac[j++] = mu * PI;
		//tinh dao ham rieng cho a7
		jac[j++] = mu * (LL - WW);
		//tinh dao ham rieng cho a8
		jac[j++] = mu * (WOPT - WW);
		//tinh dao ham rieng cho a9
		jac[j++] = mu * (N200 - A10);
		//tinh dao ham rieng cho a10
		jac[j++] = mu * (-A9);
		//tinh dao ham rieng cho b11
		jac[j++] = x * pow(XM3, B2);
		//tinh dao ham rieng cho b12
		jac[j++] = x * (WW - WOPT) * pow(XM3, B2);
		//tinh dao ham rieng cho b2
		jac[j++] = x * (B11 + B12 * (WW - WOPT)) * pow(XM3, B2) * log(XM3);
		//tinh dao ham rieng cho b3
		jac[j++] = x * pow((SS / 100.0), B4);
		//tinh dao ham rieng cho b4
		jac[j++] = x * B3 * pow((SS / 100.0), B4) * log(SS / 100.0);
		//tinh dao ham rieng cho b5
		jac[j++] = x * pow(QU, B6);
		//tinh dao ham rieng cho b6
		jac[j++] = x * B5 * pow(QU, B6) * log(QU);
		//tinh dao ham rieng cho b7
		jac[j++] = x * PI;
		//tinh dao ham rieng cho b8
		jac[j++] = x * LL;
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
			g[j]+=(f[i]-MR/PA)*jac[i][j];
	}

}
double TinhF(double* p, int numparam, int numrow)
{
	double f[NUMROW], F = 0;
	Tinhf(p,f,numparam,numrow);
	int i;
	for(i =0;i<numrow;i++)
	{
		F+= sqr(f[i]-MR/PA);
	}
	return F/2.0;
}
void xuatm(double g[], int m)
{
	for(int j=0;j<m;j++)
		printf("_a[%d]=  %.7g\n", j,g[j]);
	printf("\n");
}
#define ALPHA  0.0001
#define BETA1  0.9 
#define BETA2  0.99

void Adam(double * p, int numparam, int numrow, int maxiter=500000,
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
            v[j] = beta2t * v[j] + (1 - beta2t) * sqr(g[j]);
            mhat[j] = m[j] / (1 - beta1t);
            vhat[j] = v[j] / (1 - beta2t);
		}
		beta1t = pow(beta1, i);
        beta2t = pow(beta2, i);
        alphat = alpha * sqrt(1 - beta2t) / (1 - beta1t);
		for(j=0;j<numparam;j++)
			p[j] = p[j] - alphat * mhat[j] / (sqrt(vhat[j])+ eps);
		i++;
		//xuatm(g,numparam);
		F = TinhF(p,numparam,numrow);
		printf("i = %d => F =  %.7g\n",i, F);
	}
	xuatm(p,numparam);

}

void main()
{
	const char * filename = "F:\\_BAIBAOLEMA_Adam\\kq4.txt";
	if (!ReadData(filename, NUMROW, NUMCOL, gdata))
		printf("Can not read file %s\n", filename), exit(-1);

	double p[NUMPARAM] = { 2,13,1,-100,19,
		1,7,4,17,-0.05,3,0.000007,0.000002,2.5,
		0.3,18,4.2,-0.712,-0.004,-0.001 }; 

	Adam(p,NUMPARAM, NUMROW);
	xuatm(p,NUMPARAM);

}