#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "emissionfunctions.h"

double iZeta(double zetaPrev, double epsi, double R, double zetaSum)	{
	double zeta = zetaPrev*pow(uniform(0,1),1./R);	
	if (zeta > epsi)	{ //TODO: this was n > epsi for some reason? double check this is ok
		return iZeta(zeta, epsi, R, zetaSum+zeta);	
	}
	else	{
		return zetaSum;
	}
}

void multiGluonMC(double epsi, double R, unsigned long long N, double *I, double *stddev)	{
	*I = 0;
	double currentProductSum, r, I2 = 0;
	unsigned long long i;
	unsigned n;
	srand(time(0)); //seed RNG

	for (i=0; i < N; ++i)	{
		r = iZeta(1, epsi, R, 0);
		*I+=pow(1+r, -R)/N;
		I2+=pow(1+r, -2.*R)/N;
	}

	*stddev = pow((1./N)*(I2-pow(*I,2)), 0.5);
}

//int main(int argc, char **argv)	{
//	int i;
//	size_t x = 20;
//	unsigned long long N = 1000000;
//	double epsi = 1e-5;
//	double *R = arange(0.1,2,x);
//	
//	for(i=0; i <= x; i++)	{
//		printf("R: %f, I: %f, std: %f\n",R[i],multiGluonMC(epsi, R[i], N),stddev);	
//	}
//
//	return 0;
//}
