#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double *arange(double a, double b, size_t x)	{
	double *arr = malloc(x*sizeof(double));	
	double step = (b-a)/x;
	int i;
	
	for (i=0; i <= x; i++)	{
		arr[i] = a+step*i;
	}

	return arr;
}

double uniform(double a, double b)  {
   double range = b-a;
   return a + (rand()*range)/RAND_MAX;
}

unsigned iZeta(double zetaPrev, double epsi, double R, double *res, unsigned n)	{
	double zeta = zetaPrev*pow(uniform(0,1),1./R);	
	if (zeta > epsi)	{
		*res += zeta;
		return iZeta(zeta, epsi, R, res, n+1);	
	}
	else	{
		return n;
	}
}

double stddev;

double monteCarlo(double epsi, double R, unsigned long long N)	{
	double currentProductSum, r, I = 0, I2 = 0;
	double res;
	unsigned long long i;
	unsigned n;
	srand(time(0)); //seed RNG

	for (i=0; i < N; ++i)	{
		res = 0; //holds sum of zeta_i
		n = iZeta(1, epsi, R, &res, 0);
		I+=pow(1+res, -R)/N;
		I2+=pow(1+res, -2*R)/N;
	}

	double mu = I/N;
	stddev = pow((1./N)*(I2-pow(I,2)), 0.5);

	return I;
}

int main(int argc, char **argv)	{
	int i;
	size_t x = 20;
	unsigned long long N = 1000000;
	double epsi = 1e-5;
	double *R = arange(0.1,2,x);
	
	for(i=0; i <= x; i++)	{
		printf("R: %f, I: %f, std: %f\n",R[i],monteCarlo(epsi, R[i], N),stddev);	
	}

	return 0;
}
