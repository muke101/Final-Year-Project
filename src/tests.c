#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "emissionfunctions.h"
#include "tests.h"

double *arange(double a, double b, int x)    {
    double *arr = malloc((size_t)(x+1)*sizeof(double));
    double step = (b-a)/x;
    int i;

    for (i=0; i <= x; i++)  {
        arr[i] = a+step*i;
    }

    return arr;
}

void caMC(double x, double epsi, double R, unsigned long long N, double *I) {

    unsigned long long i;
    double phi_k, k, t, zetaV, fcorrel, zetaSum, fc, fcZ1, fcZ2, Cab, r;
	double u, u2, phi, z, z1, z2;
	struct vars transedVars;
    *I = 0;

    for (i=0; i < N; ++i)   {
        zetaV = uniform(epsi,1-epsi);
        phi_k = uniform(epsi,1-epsi);
        k = uniform(epsi,1-epsi);
        t = uniform(epsi,1-epsi);
		transedVars = transform(t, phi_k, k);
		u = transedVars.u;
		u2 = transedVars.u2;
		phi = transedVars.phi;
		z = transedVars.z;
		z1 = transedVars.z1;
		z2 = transedVars.z2;
        zetaSum = iZeta(ZETA_0, epsi, R, 0);
		fc = fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z);
		fcZ1 = fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z1);
		fcZ2 = fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z2);
        Cab = caEquation(x,k,t,u,u2,phi,z,z1,z2,fc,fcZ1,fcZ2);
        *I += Cab/N;
    }


}


void nfMC(double x, double epsi, double R, unsigned long long N, double *I) {

	*I=0;
    unsigned long long i;
    double phi_k, k, t, zetaV, fcorrel, zetaSum, fc, Cab, r, I2=0;
	double u,u2,phi,z;
	struct vars transedVars;

    for (i=0; i < N; ++i)   {
        zetaV = uniform(epsi,1-epsi);
        phi_k = uniform(epsi,1-epsi);
        k = uniform(epsi,1-epsi);
        t = uniform(epsi,1-epsi);
		transedVars = transform(t, phi_k, k);
		u = transedVars.u;
		u2 = transedVars.u2;
		phi = transedVars.phi;
		z = transedVars.z;
        zetaSum = iZeta(ZETA_0, epsi, R, 0);
		fc = fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z);
        Cab = nfEquation(x,k,u,u2,phi,z,fc); 
        *I += Cab/N;
    }


}

void caComp(double x, double epsi, double R, unsigned long long N, double *I)   {
    *I = 0;
    unsigned long long i;
    double phi_k, k, t, zetaSum, zetaV, fc, fcZ1, fcZ2, r, ca;
	double u,u2,phi,z,z1,z2;
	struct vars transedVars;

    for (i=0; i < N; ++i)   {
		phi_k = uniform(epsi,1-epsi);
		k = uniform(epsi,1-epsi);
		t = uniform(epsi,1-epsi);
		transedVars = transform(t, phi_k, k);
		u = transedVars.u;
		u2 = transedVars.u2;
		phi = transedVars.phi;
		z = transedVars.z;
		z1 = transedVars.z1;
		z2 = transedVars.z2;
        zetaSum = iZeta(ZETA_0, epsi, R, 0);
		fc = fcCorrection(z,phi,u,u2,x);
		fcZ1 = fcCorrection(z1,phi,u,u2,x);
		fcZ2 = fcCorrection(z2,phi,u,u2,x);
        ca = caEquation(x,k,t,u,u2,phi,z,z1,z2,fc,fcZ1,fcZ2);
        r = -pow(1+zetaSum, -R)*((1-x)*(11./12.*pow(M_PI,2)/6) + ca);
        *I+=r/N;
    }

}


void nfComp(double x, double epsi, double R, unsigned long long N, double *I)   {
    *I = 0;
    unsigned long long i;
    double phi_k, k, t, zetaSum, zetaV, fc, r, nf;
	double u,u2,phi,z;
	struct vars transedVars;

    for (i=0; i < N; ++i)   {
		phi_k = uniform(epsi,1-epsi);
		k = uniform(epsi,1-epsi);
		t = uniform(epsi,1-epsi);
		transedVars = transform(t, phi_k, k);
		u = transedVars.u;
		u2 = transedVars.u2;
		phi = transedVars.phi;
		z = transedVars.z;
        zetaSum = iZeta(ZETA_0, epsi, R, 0);
		fc = fcCorrection(z,phi,u,u2,x);
        nf = nfEquation(x,k,u,u2,phi,z,fc);
        r = -pow(1+zetaSum,-R)*((1-x)*(-2./12.*pow(M_PI,2)/6) + nf);
        *I+=r/N;
    }


}


void test(double defX, double defEpsi, double defR, unsigned long long N, const char flags, int step, void(*MCfunc)(double, double, double, unsigned long long, double*), void(*CompFunc)(double, double, double, unsigned long long, double*))	{
	int i,j,k;
	double I1, I2;
	double *x, *epsi, *R;
	srand(time(0));

	switch (flags)	{
		case TEST_X:
			x = arange(0,2,step);
			epsi = &defEpsi;
			R = &defR;
			for (i=0; i <= step; i++)	{
				MCfunc(x[i],*epsi,*R,N,&I1);
				CompFunc(x[i],*epsi,*R,N,&I2);
				if (x[i] == 1)	{
					printf("x=1: ");
				}
				printf("x: %f, MC: %f, Comp: %f\n", x[i], I1, I2);
			}
			free(x);
			break;
		case TEST_R:
			x = &defX;
			epsi = &defEpsi;
			R = arange(0.1,2,step);
			for (i=0; i <= step; i++)	{
				MCfunc(*x,*epsi,R[i],N,&I1);
				CompFunc(*x,*epsi,R[i],N,&I2);
				printf("R: %f, MC: %f, Comp: %f\n", R[i], I1, I2);
			}
			free(R);
			break;
		case TEST_EPSI:
			x = &defX;
            epsi = arange(1e-10,1e-3,step);
            R = &defR;
            for (i=0; i <= step; i++)	{
            	MCfunc(*x,epsi[i],*R,N,&I1);
            	CompFunc(*x,epsi[i],*R,N,&I2);
            	printf("epsi: %f, MC: %f, Comp: %f\n", epsi[i], I1, I2);
            }
			free(epsi);
			break;
		case TEST_XR:
			x = arange(0,2,step);
			epsi = &defEpsi;	
			R = arange(0.1,2,step);
			for (i=0; i <= step; i++)	{
				for (j=0; j <= step; j++)	{
					MCfunc(x[i],*epsi,R[j],N,&I1);
					CompFunc(x[i],*epsi,R[j],N,&I2);
					printf("x: %f, R: %f, MC: %f, Comp: %f\n",x[i],R[j],I1,I2);
				}
			}
			free(x);
			free(R);
			break;
		case TEST_XEPSI:
			x = arange(0,2,step);
			epsi = arange(1e-10,1e-3,step);
			R = &defR; 
			for (i=0; i <= step; i++)	{
				for (j=0; j <= step; j++)	{
					MCfunc(x[i],epsi[j],*R,N,&I1);
					CompFunc(x[i],epsi[j],*R,N,&I2);
					printf("x: %f, epsi: %f, MC: %f, Comp: %f\n",x[i],epsi[j],I1,I2);
				}
			}
			free(x);
			free(epsi);
			break;
		case TEST_REPSI:
			x = &defX;
			epsi = arange(1e-10,1e-3,step);
			R = arange(0.1,2,step);
			for (i=0; i <= step; i++)	{
				for (j=0; j <= step; j++)	{
					MCfunc(*x,epsi[i],R[j],N,&I1);
					CompFunc(*x,epsi[i],R[j],N,&I2);
					printf("epsi: %f, R: %f, MC: %f, Comp: %f\n",epsi[i],R[j],I1,I2);
				}
			}
			free(epsi);
			free(R);
			break;
		case TEST_XREPSI:
			x = arange(0,2,step);
			epsi = arange(1e-10,1e-3,step);
			R = arange(0.1,2,step);
			for (i=0; i <= step; i++)	{
				for (j=0; j <= step; j++)	{
					for (k=0; k <= step; k++)	{
						MCfunc(x[i],epsi[j],R[k],N,&I1);
						CompFunc(x[i],epsi[j],R[k],N,&I2);
						printf("x: %f, epsi: %f, R: %f, MC: %f, Comp: %f\n",x[i],epsi[j],R[k],I1,I2);
					}
				}
			}
			free(x);
			free(epsi);
			free(R);
			break;
	}

}
