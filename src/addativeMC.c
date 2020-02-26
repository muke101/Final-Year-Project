#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "emissionfunctions.h"
#define ZETA_0 1

static double fc(double t, double phi_k, double k, double x)	{
	double u2 = k/(1-k);
	double u = pow(u2, 0.5);
	double phi = phi_k*2.*M_PI; //TODO does this require the same scaling?
	double z = t;
	
	return z*pow(ua2(z,phi,u,u2),1-x/2.)+(1-z)*pow(ub2(z,phi,u,u2),-x/2.);
}

double FcorrelMin(double fcorrel, double Cab, double zetaV, double zetaSum, double R)	{
	return 1./zetaV*Cab*zetaSum*exp(-R*log(zetaV*fcorrel+1+zetaSum))*exp(-R*log(zetaV+1+zetaSum));
}

double FcorrelMaj(double fcorrel, double Cab, double zetaSum, double R)	{
	return Cab*zetaSum*exp(-R*log(fcorrel+zetaSum))*exp(-R*log(1+zetaSum));
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
		Cab = caEquation(x,k,t,phi_k)*nfEquation(x,k,t,phi_k);
		fcorrel = fc(t,phi_k,k,x);
		r = FcorrelMin(fcorrel, Cab, zetaV, zetaSum, R) + FcorrelMaj(fcorrel, Cab, zetaSum, R); 
		*I += r; 
		I2 += pow(r,2);
	}

	*stddev = pow((1./N)*(I2-pow(*I,2)),0.5);
}


int main()	{
	srand(time(0));
	double I, stddev;
	totalMC(1,1e-5,0.2,100000,&I,&stddev);
	printf("I: %f, sttdev: %f\n", I, stddev);
	return 0;
}
