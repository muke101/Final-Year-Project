#include "emissionfunctions.h"

double iZeta(double zetaPrev, double epsi, double R, double zetaSum)	{
	double zeta = zetaPrev*pow(uniform(0,1),1./R);	
	if (zeta > epsi)	{ 
		return iZeta(zeta, epsi, R, zetaSum+zeta);
	}
	else	{
		return zetaSum;
	}
}

void multiGluonMC(double epsi, double R, unsigned long long N, double *I, double *stddev)	{
	*I = 0;
	double r, I2 = 0;
	unsigned long long i;
	srand(time(0)); //seed RNG

	for (i=0; i < N; ++i)	{
		r = iZeta(1, epsi, R, 0);
		*I+=pow(1+r, -R)/N;
		I2+=pow(1+r, -2.*R)/N;
	}

	*stddev = pow((1./N)*(I2-pow(*I,2)), 0.5);
}
