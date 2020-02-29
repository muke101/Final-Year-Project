#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "emissionfunctions.h"
#define ZETA_0 1

static double fc(double t, double phi_k, double k, double x)	{
	double u2 = k/(1-k);
	double u = pow(u2, 0.5);
	double phi = phi_k*2.*M_PI; 
	double z = t;
	
	return pow(1+u2,(x-1)/2)*z*pow(ua2(z,phi,u,u2),1-x/2.)+(1-z)*pow(ub2(z,phi,u,u2),-x/2.);
}

double FcorrelMin(double fcorrel, double Cab, double zetaV, double zetaSum, double R)	{
	return 1./zetaV*Cab*(exp(-R*log(zetaV*fcorrel+1+zetaSum))-exp(-R*log(zetaV+1+zetaSum)));
}

double FcorrelMaj(double fcorrel, double Cab, double zetaSum, double R)	{
	return Cab*(exp(-R*log(fcorrel+zetaSum))-exp(-R*log(1+zetaSum)));
}

void totalMC(double x, double epsi, double R, unsigned long long N, double *I, double *stddev)	{

	unsigned long long i;
	double phi_k, k, t, zetaV, fcorrel, zetaSum, Cab, r, I2=0;
	*I = 0;
	
	for (i=0; i < N; ++i)	{
		while ((zetaV = uniform(0,1)) == 0 || zetaV == 1);
		while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
		while ((k = uniform(0,1)) == 0 || k == 1);
		while ((t = uniform(0,1)) == 0 || t == 1);
		zetaSum = iZeta(ZETA_0, epsi, R, 0);
		Cab = caEquation(x,k,t,phi_k)+nfEquation(x,k,t,phi_k);
		fcorrel = fc(t,phi_k,k,x);
		r = FcorrelMin(fcorrel, Cab, zetaV, zetaSum, R) + FcorrelMaj(fcorrel, Cab, zetaSum, R); 
		*I += r/N; 
		I2 += pow(r,2)/N;
	}

	*stddev = pow((1./N)*(I2-pow(*I,2)),0.5);
}

void caMC(double x, double epsi, double R, unsigned long long N, double *I)	{

	unsigned long long i;
	double phi_k, k, t, zetaV, fcorrel, zetaSum, Cab, r; 
	*I = 0;
	
	for (i=0; i < N; ++i)	{
		while ((zetaV = uniform(0,1)) == 0 || zetaV == 1);
		while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
		while ((k = uniform(0,1)) == 0 || k == 1);
		while ((t = uniform(0,1)) == 0 || t == 1);
		zetaSum = iZeta(ZETA_0, epsi, R, 0);
		Cab = caEquation(x,k,t,phi_k);
		fcorrel = fc(t,phi_k,k,x);
		r = FcorrelMin(fcorrel, Cab, zetaV, zetaSum, R) + FcorrelMaj(fcorrel, Cab, zetaSum, R); 
		*I += r/N; 
	}

}

void nfMC(double x, double epsi, double R, unsigned long long N, double *I)	{

	unsigned long long i;
	double phi_k, k, t, zetaV, fcorrel, zetaSum, Cab, r, I2=0;
	*I = 0;
	
	for (i=0; i < N; ++i)	{
		while ((zetaV = uniform(0,1)) == 0 || zetaV == 1);
		while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
		while ((k = uniform(0,1)) == 0 || k == 1);
		while ((t = uniform(0,1)) == 0 || t == 1);
		zetaSum = iZeta(ZETA_0, epsi, R, 0);
		Cab = nfEquation(x,k,t,phi_k);
		fcorrel = fc(t,phi_k,k,x);
		r = FcorrelMin(fcorrel, Cab, zetaV, zetaSum, R) + FcorrelMaj(fcorrel, Cab, zetaSum, R); 
		*I += r/N; 
	}

}

void caComp(double x, double epsi, double R, unsigned long long N, double *I)	{
	*I = 0;
	unsigned long long i;
	double phi_k, k, t, zetaSum, r, ca; 

	for (i=0; i < N; ++i)	{
		while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
		while ((k = uniform(0,1)) == 0 || k == 1);
		while ((t = uniform(0,1)) == 0 || t == 1);
		zetaSum = iZeta(ZETA_0, epsi, R, 0);
		ca = caEquation(x,k,t,phi_k); 
		r = -zetaSum*((1-x)*(11./12.*pow(M_PI,2)/6) + ca);
		*I+=r/N;
	}

}	

void nfComp(double x, double epsi, double R, unsigned long long N, double *I)	{
	*I = 0;
	unsigned long long i;
	double phi_k, k, t, zetaSum, r, nf; 

	for (i=0; i < N; ++i)	{
		while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
		while ((k = uniform(0,1)) == 0 || k == 1);
		while ((t = uniform(0,1)) == 0 || t == 1);
		zetaSum = iZeta(ZETA_0, epsi, R, 0);
		nf = nfEquation(x,k,t,phi_k); 
		r = -zetaSum*((1-x)*(2./12.*pow(M_PI,2)/6) + nf);
		*I+=r/N;
	}

}	



int main()	{
	srand(time(0));
	double I, stddev;
	int i;
	int step = 20;
	double *R = arange(0.1,2,step);
	double epsi = 1e-5;
	double *x = arange(0,2,step);
	double N = 10000;
	double I1, I2;
	
	//for (i=0; i <= step; i++)	{
	//	caMC(x,epsi,R[i],N,&I1);
	//	caComp(x,epsi,R[i],N,&I2);
	//	printf("FcorrelCa: %f, CaComp: %f\n", I1, I2); 
	//}
	
	for (i=0; i <= step; i++)	{
		caMC(x[i],epsi,R[2],N,&I1);
		caComp(x[i],epsi,R[2],N,&I2);
		if (x[i] == 1)	{
			printf("x=1: ");
		}
		printf("FcorrelCa: %f, CaComp: %f\n", I1, I2); 
	}

	return 0;
}
