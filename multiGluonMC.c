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

unsigned fact(unsigned n, unsigned t)	{
	if (n == 1)	{
		return t;	
	}
	return fact(n-1,t*n);	
}

unsigned iZeta(double zetaPrev, double epsi, double R, double *res, unsigned n)	{
	double zeta = zetaPrev*pow(uniform(epsi,1),1./R);	
	if (zeta > epsi)	{
		res[0] *= 1./zeta; //this quickly sky rockets to infinity
		res[1] += zeta;
		return iZeta(zeta, epsi, R, res, n+1);	
	}
	else	{
		return n;
	}
}

double monteCarlo(double epsi, double R, unsigned long long N)	{
	double currentProductSum, r, I = 0, I2 = 0;
	double res[2];
	unsigned long long i;
	unsigned n;
	srand(time(0)); //seed RNG

	res[0] = 1; //holds product sum of 1/zeta_i
	for (i=0; i < N; ++i)	{
		currentProductSum = res[0];
		res[1] = 0; //holds sum of zeta_i
		n = iZeta(1, epsi, R, res, 0);
		if (res[1] < 1 && n > 0)	{
			r = pow(R,n)*res[0]*(1-res[1]/N)/fact(n,1);
			I+=r;
			I2+=pow(r,2);
		}
		else	{
			res[0] = currentProductSum; //undo invalidated changes to product sum
		}
	}

	return pow(epsi,R)*I;
}

int main(int argc, char **argv)	{
	int i;
	size_t x = 20;
	unsigned long long N = 10000;
	double epsi = 10e-3;
	double *R = arange(0.1,2,x);
	
	for(i=0; i <= x; i++)	{
		printf("%f\n",monteCarlo(epsi, R[i], N));	
	}

	return 0;
}
