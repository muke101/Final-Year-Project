#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "emissionfunctions.h"
#include "tests.h"

double *arange(double a, double b, int x)    {
    double *arr = malloc((size_t)x*sizeof(double));
    double step = (b-a)/x;
    int i;

    for (i=0; i <= x; i++)  {
        arr[i] = a+step*i;
    }

    return arr;
}

void caMC(double x, double epsi, double R, unsigned long long N, double *I) {

    unsigned long long i;
    double phi_k, k, t, zetaV, fcorrel, zetaSum, Cab, r;
	double u, u2, phi, z, z1, z2;
	double *transedVars;
    *I = 0;

    for (i=0; i < N; ++i)   {
        while ((zetaV = uniform(0,1)) == 0 || zetaV == 1);
        while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
        while ((k = uniform(0,1)) == 0 || k == 1);
        while ((t = uniform(0,1)) == 0 || t == 1);
		transedVars = transform(t, phi_k, k);
		u = transedVars[0];
		u2 = transedVars[1];
		phi = transedVars[2];
		z = transedVars[3];
		z1 = transedVars[4];
		z2 = transedVars[5];
        zetaSum = iZeta(ZETA_0, epsi, R, 0);
        Cab = caEquation(x,k,t,u,u2,phi,z,z1,z2,fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z),fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z1),fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z2));
        fcorrel = fc(u,u2,phi,z,x);
        r = FcorrelMin(fcorrel, Cab, zetaV, zetaSum, R) + FcorrelMaj(fcorrel, Cab, zetaSum, R);
        *I += r/N;
    }

}


void nfMC(double x, double epsi, double R, unsigned long long N, double *I) {

	*I=0;
    unsigned long long i;
    double phi_k, k, t, zetaV, fcorrel, zetaSum, Cab, r, I2=0;
	double u,u2,phi,z;
	double *transedVars;

    for (i=0; i < N; ++i)   {
        while ((zetaV = uniform(0,1)) == 0 || zetaV == 1);
        while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
        while ((k = uniform(0,1)) == 0 || k == 1);
        while ((t = uniform(0,1)) == 0 || t == 1);
		transedVars = transform(t, phi_k, k);
		u = transedVars[0];
		u2 = transedVars[1];
		phi = transedVars[2];
		z = transedVars[3];
        zetaSum = iZeta(ZETA_0, epsi, R, 0);
        Cab = nfEquation(x,k,u,u2,phi,z,fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z)); //TODO is this fcVsc or fcCorrection?
        fcorrel = fc(u,u2,phi,z,x);
        r = FcorrelMin(fcorrel, Cab, zetaV, zetaSum, R) + FcorrelMaj(fcorrel, Cab, zetaSum, R);
        *I += r/N;
    }

}

void caComp(double x, double epsi, double R, unsigned long long N, double *I)   {
    *I = 0;
    unsigned long long i;
    double phi_k, k, t, zetaSum, zetaV, r, ca;
	double u,u2,phi,z,z1,z2;
	double *transedVars;

    for (i=0; i < N; ++i)   {
		while ((zetaV = uniform(0,1)) == 0);
        while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
        while ((k = uniform(0,1)) == 0 || k == 1);
        while ((t = uniform(0,1)) == 0 || t == 1);
		transedVars = transform(t, phi_k, k);
		u = transedVars[0];
		u2 = transedVars[1];
		phi = transedVars[2];
		z = transedVars[3];
		z1 = transedVars[4];
		z2 = transedVars[5];
        zetaSum = iZeta(ZETA_0, epsi, R, 0);
        ca = caEquation(x,k,t,u,u2,phi,z,z1,z2,fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z),fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z1),fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z2));
        r = -pow(1+zetaSum, R)*((1-x)*(11./12.*pow(M_PI,2)/6) + ca);
        *I+=r/N;
    }

}


void nfComp(double x, double epsi, double R, unsigned long long N, double *I)   {
    *I = 0;
    unsigned long long i;
    double phi_k, k, t, zetaSum, zetaV, r, nf;
	double u,u2,phi,z;
	double *transedVars;

    for (i=0; i < N; ++i)   {
		while ((zetaV = uniform(0,1)) == 0);
        while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
        while ((k = uniform(0,1)) == 0 || k == 1);
        while ((t = uniform(0,1)) == 0 || t == 1);
		transedVars = transform(t, phi_k, k);
		u = transedVars[0];
		u2 = transedVars[1];
		phi = transedVars[2];
		z = transedVars[3];
        zetaSum = iZeta(ZETA_0, epsi, R, 0);
        nf = nfEquation(x,k,u,u2,phi,z,fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z));
        r = -zetaSum*((1-x)*(2./12.*pow(M_PI,2)/6) + nf);
        *I+=r/N;
    }

}


void test(double defX, double defEpsi, double defR, unsigned long long N, const char flags, int step, void(*MCfunc)(double, double, double, unsigned long long, double*), void(*CompFunc)(double, double, double, unsigned long long, double*))	{
	int i;
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
				printf("MC: %f, Comp: %f\n", I1, I2);
			}
			break;
		case TEST_R:
			x = &defX;
			epsi = &defEpsi;
			R = arange(0.1,2,step);
			for (i=0; i <= step; i++)	{
				MCfunc(*x,*epsi,R[i],N,&I1);
				CompFunc(*x,*epsi,R[i],N,&I2);
				printf("MC: %f, Comp: %f\n", I1, I2);
			}
			break;
		case TEST_EPSI:
			x = &defX;
            epsi = arange(1e-10,1e-3,step);
            R = &defR;
            for (i=0; i <= step; i++)	{
            	MCfunc(*x,epsi[i],*R,N,&I1);
            	CompFunc(*x,epsi[i],*R,N,&I2);
            	printf("MC: %f, Comp: %f\n", I1, I2);
            }
			break;
		case TEST_XR:
			x = arange(0,2,step);
			epsi = &defEpsi;	
			R = arange(0.1,2,step);
			for (i=0; i <= step; i++)	{
				MCfunc(x[i],*epsi,R[i],N,&I1);
				CompFunc(x[i],*epsi,R[i],N,&I2);
				printf("MC: %f, Comp: %f\n", I1, I2);
			}
			break;
		case TEST_XEPSI:
			x = arange(0,2,step);
			epsi = arange(1e-10,1e-3,step);
			R = &defR; 
			for (i=0; i <= step; i++)	{
				MCfunc(x[i],epsi[i],*R,N,&I1);
				CompFunc(x[i],epsi[i],*R,N,&I2);
				printf("MC: %f, Comp: %f\n", I1, I2);
			}
			break;
		case TEST_REPSI:
			x = &defX;
			epsi = arange(1e-10,1e-3,step);
			R = arange(0.1,2,step);
			for (i=0; i <= step; i++)	{
				MCfunc(*x,epsi[i],R[i],N,&I1);
				CompFunc(*x,epsi[i],R[i],N,&I2);
				printf("MC: %f, Comp: %f\n", I1, I2);
			}
			break;
		case TEST_XREPSI:
			x = arange(0,2,step);
			epsi = arange(1e-10,1e-3,step);
			R = arange(0.1,2,step);
			for (i=0; i <= step; i++)	{
				MCfunc(x[i],epsi[i],R[i],N,&I1);
				CompFunc(x[i],epsi[i],R[i],N,&I2);
				printf("MC: %f, Comp: %f\n", I1, I2);
			}
			break;
	}

}
