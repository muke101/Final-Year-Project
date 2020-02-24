#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define ZETA_0 1


void totalMC(double x, double epsi, double R, unsigned long long N, double *I, double *stddev)	{
	unsigned long long i;
	double mu, phi_k, k, t, zetaV, fcorrel, zetaSum, I2=0;
	*I = 0;
	
	for (i=0; i < N; ++i)	{
		phi_k = uniform(0,1);
		t = uniform(0,1);
		k = uniform(0,1);
		zetaV = uniform(0,1);
		fcorrel = fc(x, k, t, phi_k);
		zetaSum = iZeta(ZETA_0, epsi, R);
		*I += FcorrelMin(fcorrel, zetaV, zetaSum) + FcorrelMaj(fcorrel, zetaSum);
		I2 += pow(I,2);
	}

	mu = I/N;
	*stddev = pow((1./N)*(I2-pow(I,2)),0.5);
}


int main()	{
	return 0;
}
